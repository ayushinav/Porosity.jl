# Fixed discretization

Geophysical models generally have fixed discretization. This is mostly because the different numerical schemes such as finite difference and finite element are computationally expensive and allocating a grid prior to solving the corresponding PDEs saves some computational resources. We provide the capability to do MCMC inference on such fixed grids.

Let's denote the model parameters, eg., conductivity, by `m`, and the layer thickness by `h`. Therefore, in a N-layer case, we will have

```math
m = [m_1, m_2, m_3, ... , m_N] \\
h = [h_1, h_2, h_3, ... , h_{N-1}]
```

such that

```math
m\_i  \in \mathcal{D}_{m_i} \text{ ; where } \mathcal{D}_{m_i} = \textit{a priori} \text{ distribution for } m_i
```

and $h_i$ is fixed.

In the following example, we show how to perform MCMC inversion for such a case using a synthetic dataset. We assume a 6-layered earth, including the half-space where all layers have resistivities bounded between $10^{-1}$ and $10^5$, defined using a uniform distribution.

## Copy-Pasteable code

```@example fixed_mcmc
using MT
using Distributions
using Turing
using LinearAlgebra
using CairoMakie

m_test = MTModel(log10.([100.0, 10.0, 1000.0]), [1e3, 1e3]);
f = 10 .^ range(-4; stop=1, length=25);
ω = vec(2π .* f);

r_obs = forward(m_test, ω);

err_phi = asin(0.02) * 180 / π .* ones(length(ω));
err_appres = 0.1 * r_obs.ρₐ;
err_resp = MTResponse(err_appres, err_phi);

r_obs.ρₐ .= r_obs.ρₐ .+ randn(size(f)) .* err_appres;
r_obs.ϕ .= r_obs.ϕ .+ randn(size(f)) .* err_phi;

respD = MTResponseDistribution(normal_dist, normal_dist);

z = collect(0:500:2.5e3)
h = diff(z);

# fixed discretization
modelD = MTModelDistribution(
    Product(
        [Uniform(-1.0, 5.0) for i in eachindex(z)]
    ),
    vec(h)
);

n_samples = 1000;
mcache = mcmc_cache(modelD, respD, n_samples, NUTS());

mt_chain = stochastic_inverse(r_obs, err_resp, ω, mcache)
```

The obtained `mt_chain` contains the *a posteriori* distributions that can be saved using [JLD2.jl](https://github.com/JuliaIO/JLD2.jl).

```julia
using JLD2
JLD2.@save "file_path.jld2" mt_chain
```

and plotted as :

```@example fixed_mcmc
fig = Figure()
ax = Axis(fig[1, 1]; xscale=log10)
hm = get_kde_image!(ax, mt_chain, modelD; kde_transformation_fn=log10,
    trans_utils=(; m=pow_tf), colormap=:thermal, colorrange=(-3.0, 0))
Colorbar(fig[1, 2], hm; label="log pdf")

mean_kws = (; color=:blue, linewidth=2)
std_kws = (; color=:red, linewidth=2)
get_mean_std_image!(ax, mt_chain, modelD; confidence_interval=0.9, mean_kwargs=mean_kws,
    std_plus_kwargs=std_kws, std_minus_kwargs=std_kws)
xlims!(ax, [1e-1, 1e5])
ylims!(ax, [2500, 0])

plot_model!(ax, m_test; color=:black, linestyle=:dash, label="true")
Legend(fig[2, :], ax; orientation=:horizontal)
fig
```

The list of models can then be obtained from chains using

```@example fixed_mcmc
model_list = get_model_list(mt_chain, modelD);
nothing # hide
```

We can then easily check the fit of the response curves

```@example fixed_mcmc
fig = Figure()
ax1 = Axis(fig[1, 1])
ax2 = Axis(fig[1, 2])

resp_post = forward(model_list[1], ω);
for i in 1:(length(model_list) > 100 ? 100 : length(model_list))
    forward!(resp_post, model_list[i], ω)
    plot_response!([ax1, ax2], ω, resp_post; alpha=0.4, color=:gray)
end

plot_response!([ax1, ax2], ω, r_obs; errs=err_resp, plt_type=:errors, whiskerwidth=10)
plot_response!([ax1, ax2], ω, r_obs; plt_type=:scatter, label="true")

fig
```

The posterior distribution can then be obtained as:

```julia
pre_img = pre_image(m_dist, mt_chain);
kde_img = get_kde_image(pre_img..., false; xscale=:identity, yscale=:identity, yflip=true)
```

We can also obtain the mean and 1 std deviation bounds as:

```julia
mean_std_plt_lin = get_mean_std_image(pre_img...; yscale=:identity)
```
