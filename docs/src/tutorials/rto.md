## RTO-TKO

RTO-TKO is a stochastic framework introduced in the electromagnetic geophysical community by Blatter et al., 2022 ([a](https://doi.org/10.1093/gji/ggac241) and [b](https://doi.org/10.1093/gji/ggac242))

Similar to fixed discretization scheme, the grid sizes do not change. RTO-TKO chooses a perturbation in the prior space and optimzes for a response sampled from the response domain.

Let's have

```math
m = [m_1, m_2, m_3, ... , m_N] \\
h = [h_1, h_2, h_3, ... , h_{N_1}]
```

and

```math
d = \mathcal{F}(m)
```

Let ``C_d`` be the data covariance matrix and we want to explore uncertainty for the following misfit function:

```math
J(m) = [\mathcal{F}(m) - d]^T C_d [\mathcal{F}(m) - d] + \mu (Lm)^T(Lm)
```

where ``\mu`` is the regularization weight and ``L`` is the derivative matrix. RTO-TKO explores the uncertainty in ``\mu``-space instead of fixing it and gives a family of models that fit the data. The algorithm works as:

Solving for

```math
\begin{align*}
& J(m) = [\mathcal{F}(m) - d]^T C_d [\mathcal{F}(m) - d] + \mu (Lm)^T(Lm) \\

& 1) \quad \text{Solve for } \; m_{i+1} \\
& \text{Sample  } \; \tilde{d} \sim \mathcal{N}(d, C_d) \text{ and } \tilde{m} \sim \mathcal{N}(0, \frac{1}{\mu}(L^T L))\\

& \text{Solve} \\
& m_{i+1} = \argmin_{m_{i+1}} \quad [\mathcal{F}(m_{i+1}) - \tilde{d}]^T C_d [\mathcal{F}(m_{i+1}) - \tilde{d}] + \mu_i [L(m_{i+1} - \tilde{m})]^T[L(m_{i+1} - \tilde{m})] \\

& 2) \quad \text{Solve for } \; \mu_{i+1} \\
& \text{Sample  } \; \tilde{d} \sim \mathcal{N}(d, C_d) \\

& \text{Solve} \\
& \mu_{i+1} = \argmin_{\mu_{i+1}} \quad [\mathcal{F}(m_{i+1}) - \tilde{d}]^T C_d [\mathcal{F}(m_{i+1}) - \tilde{d}] - log(p(\mu_{i+1})) \\
\end{align*}
```

!!! note
    
      - Usually, the prior of ``\mu``is a uniform distribution and we do not have to compute the log pdf term
      - Implementing the above from scratch might not be trivial because of ``L'L`` being non-invertible, and we do the optimization in the domain defined by $\xi = \sqrt{\mu}Lm$.
        Such a variable will then have a standard normal distribution $\mathcal{N}(0, I)$ when $m \sim \mathcal{N}(0, \frac{1}{\mu}(L^T L))$

## Copy-Pasteable code

```@example rto_tko
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

z = collect(0:100:2.5e3)
h = diff(z)

# obtaining Occam result since RTO-TKO is closely derived from it
m_occam = MTModel(2 .* ones(length(z)), vec(h))
W = diagm(inv.([err_resp.ρₐ..., err_resp.ϕ...])) .^ 2

ret_code = inverse!(
    m_occam,
    r_obs,
    ω,
    Occam(; μgrid=[1e-2, 1e6]);
    W=W,
    χ2=1.0,
    max_iters=50,
    verbose=true
)

respD = MTResponseDistribution(normal_dist, normal_dist);

modelD = MTModelDistribution(
    product_distribution([Uniform(-1.0, 5.0) for i in eachindex(z)]),
    vec(h)
)

n_samples = 100
m_rto = MTModel(2 .* ones(length(z)), vec(h));
r_cache = rto_cache(
    m_rto, [1e-2, 1e3], Occam(), n_samples, n_samples, 1.0, [:ρₐ, :ϕ], false);

rto_chain = stochastic_inverse(
    r_obs,
    err_resp,
    ω,
    r_cache;
    model_trans_utils=(; m=MT.lin_tf)
)

mt_chain = Turing.Chains(
    (rto_chain.value.data[:, 1:(end - 1), :]),
    [Symbol("m[$i]") for i in 1:length(z)]
)
```

The obtained `mt_chain` contains the distributions that can be saved using [JLD2.jl](https://github.com/JuliaIO/JLD2.jl).

```julia
using JLD2
JLD2.@save "file_path.jld2" mt_chain
```

and plotted as :

```@example rto_tko
fig = Figure()
ax = Axis(fig[1, 1]; xscale=log10)
hm = get_kde_image!(ax, mt_chain, modelD; kde_transformation_fn=log10,
    trans_utils=(; m=pow_tf), colormap=:thermal, colorrange=(-2.0, -1.5))
Colorbar(fig[1, 2], hm; label="log pdf")

mean_kws = (; color=:blue, linewidth=2)
std_kws = (; color=:red, linewidth=2)
get_mean_std_image!(ax, mt_chain, modelD; confidence_interval=0.9, mean_kwargs=mean_kws,
    std_plus_kwargs=std_kws, std_minus_kwargs=std_kws)
ylims!(ax, [2500, 0])

plot_model!(ax, m_test; color=:black, linestyle=:dash, label="true")
plot_model!(ax, m_occam; color=:green, label="occam")
Legend(fig[2, :], ax; orientation=:horizontal)
fig
```

The list of models can then be obtained from chains using

```@example rto_tko
model_list = get_model_list(mt_chain, modelD);
nothing # hide
```

We can then easily check the fit of the response curves

```@example rto_tko
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
ylims!(ax1, exp10.([1.2, 3.1]))
ylims!(ax2, [12, 70])

fig
```

!!! note
    
    Samples from RTO-TKO can be usually filtered out to reject the samples that have a poor fit on the data. While we do not filter them here for the sake of brevity, we have explicitly defined the axes limits to exclude their response curves. Some of them are still there on the plot. This is done to make user aware that such samples might exist.
