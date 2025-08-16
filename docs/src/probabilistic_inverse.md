# Probabilistic inversion

Performing stochastic inversion involves forming the model prior space, defining the likelihood and then getting samples from the posterior space.

## Constructing distributions to sample from

### Model distribution (*a priori* information)

Before beginning to talk about how to construct the *a priori* distribution, it is important to understand that the model here refers to the discretization space as well as the values of physical properties. For eg, for electrical methods, the model space will consist of the electrical resistivity as well as the grid sizes represented by those. This is primarily useful in 1D, but can have its own consequences.

For most applications, however, we fix the node points. We follow a 1D MT example to show the framework. For now, we begin by choosing a very broad prior for all the `n` layers. A model distribution can be constructed using `MTModelDistribution(...)`. This has the same structure as `MTModel`, that is, the first parameter denotes the prior for electrical conductivities, while the second is for the layer thicknesses `h`. If you do not want to infer on `h`, just pass it as a simple vector, as is also demonstrated in the following case:

```julia
n = 50 # number of layers
h = fill(20.0, n - 1) # making the discretization

# prior space with fixed model discretization
modelD = MTModelDistribution(
    Product(
        [Uniform(-1, 5) for i in 1:n]
    ),
    vec(h) # fixed h
);
```

Some of the things to understand here are :

  - `Product(...)` is a `Distributions.jl` function that joins a series of multivariate/univariate samples, here `Uniform(-1,5)`. Thus, `Product([Uniform(-1, 5) for in 1:n])` makes a sampler that will output an `n`-length vector, with each element sampled from `Uniform(-1,5)`. Here, we use a uniform prior, but one can choose any prior. The only thing to keep in mind is that the prior can be sampled, that is, `rand(modelD.m)` returns a vector of appropriate length. Another eg., would be `MultivariateNormal(μ, Σ)`.
  - To also sample `h`, pass another multivariate distribution that outputs `n-1` length vector. Again using a uniform, independent prior, we'll have:

```julia

# prior space with variable model discretization
modelD = MTModelDistribution(
    Product(
        [Uniform(-1, 5) for i in 1:n]
    ),
    Product(
        [Uniform(15, 25) for i in 1:(n - 1)]
    )
);
```

This has its own consequences on the results and one needs to be mindful of that while interpreting the results.

### Response distribution (*likelihood*)

A likelihood is determined by an observed response `r_obs` we want to fit, and the errors associated with it `err_resp`. One of the popular ways in which likelihood can be formed is using the gaussian distribution, centered around `r_obs` with variance given by `err_resp`. To make things consistent, we pass `r_obs` and `err_resp` into the final function. To make a likelihood, we just need to pass a function that can take in a response parameter and the associated error and produce a distribution. We already provide a function `norm_dist` that takes in a vector for the response and another vector/matrix for the covariance matrix of the errors. The response distribution is then constructed by:

```julia
respD = MTResponseDistribution(normal_dist, normal_dist)
```

!!! note
    
    Both `r_obs` and `err_resp` have the same type, eg. `MTResponse`.

## Inference

We now have all the ingredients to perform inversion, except the sampler. This [page](https://turing.ml/dev/docs/using-turing/sampler-viz) provides a brief review of samplers `Turing.jl` provides. The number of posterior points to be sampled `n_samples`, the algorithm `mcmc_alg` and the distributions are brought together by `mcmc_cache`.

```julia
mcmc_alg = NUTS();
n_samples = 40;
mcache = mcmc_cache(modelD, respD, n_samples, mcmc_alg);
```

The posterior samples are then sampled by simply calling:

```julia
mcmc_chain = stochastic_inverse(r_obs, err_resp, ω, mcache)
```

## Copy-Pasteable code

```julia
using MT
using Distributions
using Turing
using LinearAlgebra

m_test = MTModel(log10.([100.0, 10.0, 1000.0]), [1e3, 1e3]);
f = 10 .^ range(-4; stop=1, length=25);
ω = vec(2π .* f);

r_obs = forward(m_test, ω);

err_phi = asin(0.01) * 180 / π .* ones(length(ω));
err_appres = 0.02 * r_obs.ρₐ;
err_resp = MTResponse(err_appres, err_phi);

r_obs.ρₐ .= r_obs.ρₐ .+ err_appres;
r_obs.ϕ .= r_obs.ϕ .+ err_phi;

respD = MTResponseDistribution(normal_dist, normal_dist);

z = 10 .^ collect(range(1; stop=4, length=100));
h = diff(z);

modelD = MTModelDistribution(
    Product(
        [Uniform(-1.0, 5.0) for i in eachindex(z)]
    ),
    vec(h)
);

n_samples = 50;
mcache = mcmc_cache(modelD, respD, 50, NUTS());

mcmc_chain = stochastic_inverse(r_obs, err_resp, ω, mcache)
```

The obtained `mcmc_chain` contains the distributions that can be saved using [JLD2.jl](https://github.com/JuliaIO/JLD2.jl).

```julia
using JLD2
JLD2.@save "file_path.jld2" mcmc_chains
```

**Note**:

!!! note
    
    The returned chains will be sampled in the distribution specified by `modelD`. In the presented case, it will have values $\in [-1, 5]$ and we can get the values by `10. ^ value`.

The list of models can then be obtained from chains using

```
model_list = get_model_list(mcmc_chains, modelD)
```
