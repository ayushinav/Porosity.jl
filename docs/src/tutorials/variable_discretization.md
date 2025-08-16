## Variable discretization

Especially in 1D, we have an analytical solution to the model with any grid spacing. This allows us to move the layer interfaces up and down while parameterizing the model space in a different way wherein there is more flexibility.

Let's denote the model parameters, eg., conductivity, by `m`, and the layer thickness by `h`. Therefore, in a N-layer case, we will have

```math
m = [m_1, m_2, m_3, ... , m_N] \\
h = [h_1, h_2, h_3, ... , h_{N_1}]
```

such that

```math
m\_i  \in \mathcal{D}_{m_i} \text{ ; where } \mathcal{D}_{m_i} = \textit{a priori} \text{ distribution for } m_i \text{ ; and} \\
 h\_i  \in \mathcal{D}_{h_i} \text{ ; where } \mathcal{D}_{h_i} = \textit{a priori} \text{ distribution for } h_i 
```

In the following example, we show how to perform MCMC inversion for such a case using a synthetic dataset. We assume a 6-layered earth, including the half-space where all layers have resistivities bounded between $10^{-1}$ and $10^5$, defined using a uniform distribution, and layered thickness values vary between 100 m and 500 m for all the layers above the half-space.

## Copy-Pasteable code

```@example variable_mcmc
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
h_ub = h .+ 100;
h_lb = h .- 100;

# fixed discretization
modelD = MTModelDistribution(
    Product(
        [Uniform(-1.0, 5.0) for i in eachindex(z)]
    ),
    Product(
        [Uniform(h_lb[i], h_ub[i]) for i in eachindex(h)]
    )
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

```@example variable_mcmc
fig = Figure()
ax = Axis(fig[1, 1]; xscale=log10)
hm = get_kde_image!(
    ax, mt_chain, modelD; kde_transformation_fn=log10, trans_utils=(; m=pow_tf),
    grid=(m=collect(-1:0.1:5), z=collect(1:50:2.5e3)),
    colormap=:thermal, colorrange=(-4, 0.0))
Colorbar(fig[1, 2], hm; label="log pdf")

mean_kws = (; color=:blue, linewidth=2)
std_kws = (; color=:red, linewidth=2)
get_mean_std_image!(
    ax, mt_chain, modelD; confidence_interval=0.9, trans_utils=(; m=pow_tf),
    mean_kwargs=mean_kws, std_plus_kwargs=std_kws,
    std_minus_kwargs=std_kws, z_points=collect(1:10:2.5e3))
xlims!(ax, [1e-1, 1e5])
ylims!(ax, [2500, 0])

plot_model!(ax, m_test; color=:black, linestyle=:dash, linewidth=2, label="true")
Legend(fig[2, :], ax; orientation=:horizontal)
fig
```

The list of models can then be obtained from chains using

```@example variable_mcmc
model_list = get_model_list(mt_chain, modelD);
nothing # hide
```

We can then easily check the fit of the response curves

```@example variable_mcmc
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
