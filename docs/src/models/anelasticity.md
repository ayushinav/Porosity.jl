# Anelastic models

```@setup anelastic_plts
using Porosity, CairoMakie, UnPack
```

## Andrade pseudoperiod model

```@docs; canonical = false
andrade_psp
```

For different temperatures, the distribution with oscillation period looks like (compare with Fig. 1 (e) and (f) of [Jackson and Faul, 2010](https://doi.org/10.1016/j.pepi.2010.09.005))

!!! tip
    
    The following code is a nice beginner example on [changing `params`](../intro/getting_started.md#Changing-initial-parameters) values for rock physics models

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example anelastic_plts
T = (900:50:1200) .+ 273.0
P = 0.2
dg = 3.1
σ = 10.0 * 1.0f-3
ϕ = 0.0
ρ = 3300.0

f = 10 .^ -collect(range(0, 3; length=100))

params = default_params(andrade_psp)
params_elastic = default_params(anharmonic)
@unpack T_K_ref, dG_dT, dG_dP, P_Pa_ref = params_elastic
new_Gu_ol = 1e-9 * (params.G_UR * 1e9 - (900 + 273 - T_K_ref) * dG_dT -
             (0.2 * 1e9 - P_Pa_ref) * dG_dP)

new_params_elastic = (; params_elastic..., Gu_0_ol=new_Gu_ol)

new_params = (; params..., params_elastic=new_params_elastic)

m = andrade_psp(T, P, dg, σ, ϕ, ρ, f')

resp = forward(m, [], new_params);
size(resp.J1)

fig = Figure()
ax = Axis(fig[1, 1]; xscale=log10, backgroundcolor=(:magenta, 0.05),
    ylabel="Shear Modulus (GPa)", xlabel="T(s)")

for i in eachindex(T)
    color_ = RGBf(i / 10, 0.0, 1 - i / 10)
    lines!(ax, 1 ./ f, resp.M[i, :] ./ (1e9); color=color_, label="$(T[i] - 273)")
end
xlims!(ax, 1e-2, 1e4)
ylims!(ax, 0, 80)

ax2 = Axis(fig[1, 2]; xscale=log10, yscale=log10,
    backgroundcolor=(:magenta, 0.05), ylabel="Q⁻¹ ", xlabel="T(s)")

for i in eachindex(T)
    color_ = RGBf(i / 10, 0.0, 1 - i / 10)
    lines!(ax2, 1 ./ f, resp.Qinv[i, :]; color=color_, label="$(T[i] - 273)")
end
xlims!(ax2, 1e-2, 1e4)
ylims!(ax2, 10.0^(-2.5), 10.0^0.5)
fig[2, 1:2] = Legend(fig, ax, "Temperature"; orientation=:horizontal, labelsize=12)

nothing # hide
```

```@raw html
</details>
```

```@example anelastic_plts
fig # hide
```

## Extended Burgers model

```@docs; canonical = false
eburgers_psp
```

For different temperatures, the distribution with oscillation period looks like (compare with Fig. 1 (a) and (b) of [Jackson and Faul, 2010](https://doi.org/10.1016/j.pepi.2010.09.005))

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example anelastic_plts
T = (600:50:1200) .+ 273.0
P = 0.2
dg = 3.1
σ = 10.0 * 1.0f-3
ϕ = 0.0
ρ = 3300.0

f = 10 .^ -collect(range(-2, 4; length=100))

m = eburgers_psp(T, P, dg, σ, ϕ, ρ, f')
resp = forward(m, []);
size(resp.J1)

fig = Figure()
ax = Axis(fig[1, 1]; xscale=log10, backgroundcolor=(:magenta, 0.05),
    ylabel="Shear Modulus (GPa)", xlabel="T(s)")

for i in eachindex(T)
    color_ = RGBf(i / 10, 0.0, 1 - i / 10)
    lines!(ax, 1 ./ f, resp.M[i, :] ./ (1e9); color=color_, label="$(T[i] - 273)")
end
xlims!(ax, 1e-2, 1e4)
ylims!(ax, 0, 80)

ax2 = Axis(fig[1, 2]; xscale=log10, yscale=log10,
    backgroundcolor=(:magenta, 0.05), ylabel="Q⁻¹ ", xlabel="T(s)")

for i in eachindex(T)
    color_ = RGBf(i / 10, 0.0, 1 - i / 10)
    lines!(ax2, 1 ./ f, resp.Qinv[i, :]; color=color_, label="$(T[i] - 273)")
end
xlims!(ax2, 1e-2, 1e4)
ylims!(ax2, 10.0^(-2.5), 10.0^0.5)
fig[2, 1:2] = Legend(
    fig, ax, "Temperature"; orientation=:horizontal, labelsize=12, nbanks=2)

nothing # hide
```

```@raw html
</details>
```

```@example anelastic_plts
fig # hide
```

!!! warn
    
    Use of `peak` coefficients is highly discouraged because these involve integrals that can be unstable in their limits. Note that the figures are skewed and the weird behaviour increases at lower temperatures.

Similar figures for different `params` can be obtained (compare with Fig. 1 (c) and (d) of [Jackson and Faul, 2010](https://doi.org/10.1016/j.pepi.2010.09.005)):

!!! tip
    
    The following code is a nice example on [changing `params`](../intro/getting_started.md#Changing-initial-parameters) values for rock physics models

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example anelastic_plts
T = (600:50:1200) .+ 273.0
P = 0.2
dg = 3.1
σ = 10.0 * 1.0f-3
ϕ = 0.0
ρ = 3300.0

params = default_params(eburgers_psp)
params_elastic = default_params(anharmonic)
@unpack T_K_ref, dG_dT, dG_dP, P_Pa_ref = params_elastic
new_Gu_ol = 1e-9 *
            (Porosity.params_JF10.s6585_bg_peak.G_UR * 1e9 - (900 + 273 - T_K_ref) * dG_dT -
             (0.2 * 1e9 - P_Pa_ref) * dG_dP)

new_params_elastic = (; params_elastic..., Gu_0_ol=new_Gu_ol)

new_params = (; params..., params_elastic=new_params_elastic,
    params_btype=Porosity.params_JF10.s6585_bg_peak)

m = eburgers_psp(T, P, dg, σ, ϕ, ρ, f')

resp = forward(m, [], new_params);
size(resp.J1)

fig = Figure(; size=(700, 400))
ax = Axis(fig[1, 1]; xscale=log10, backgroundcolor=(:magenta, 0.05),
    ylabel="Shear Modulus (GPa)", xlabel="T(s)")

for i in eachindex(T)
    color_ = RGBf(i / 10, 0.0, 1 - i / 10)
    lines!(ax, 1 ./ f, resp.M[i, :] ./ (1e9); color=color_, label="$(T[i] - 273)")
end
xlims!(ax, 1e-2, 1e4)
ylims!(ax, 0, 80)

ax2 = Axis(fig[1, 2]; xscale=log10, yscale=log10,
    backgroundcolor=(:magenta, 0.05), ylabel="Q⁻¹ ", xlabel="T(s)")

for i in eachindex(T)
    color_ = RGBf(i / 10, 0.0, 1 - i / 10)
    lines!(ax2, 1 ./ f, resp.Qinv[i, :]; color=color_, label="$(T[i] - 273)")
end
xlims!(ax2, 1e-2, 1e4)
ylims!(ax2, 10.0^(-2.5), 10.0^0.5)
fig[2, 1:2] = Legend(
    fig, ax, "Temperature"; orientation=:horizontal, labelsize=12, nbanks=2)

nothing # hide
```

```@raw html
</details>
```

```@example anelastic_plts
fig # hide
```

## Premelt model

```@docs; canonical = false
premelt_anelastic
```

For different temperatures, the distribution with oscillation period looks like

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example anelastic_plts
T = (600:50:1400) .+ 273.0
P = 0.2
dg = 3.1
σ = 10.0 * 1.0f-3
ϕ = 0.0
ρ = 3300.0
T_solidus = 1100 + 273.0

f = 10 .^ -collect(range(-2, 4; length=100))

m = premelt_anelastic(T, P, dg, σ, ϕ, ρ, 0.0, T_solidus, f')

resp = forward(m, []);
size(resp.J1)

fig = Figure()
ax = Axis(fig[1, 1]; xscale=log10, backgroundcolor=(:magenta, 0.05),
    ylabel="Shear Modulus (GPa)", xlabel="T(s)")

for i in eachindex(T)
    color_ = RGBf(i / 10, 0.0, 1 - i / 10)
    lines!(ax, 1 ./ f, resp.M[i, :] ./ (1e9); color=color_, label="$(T[i] - 273)")
end
xlims!(ax, 1e-2, 1e4)
ylims!(ax, 0, 80)

ax2 = Axis(fig[1, 2]; xscale=log10, yscale=log10,
    backgroundcolor=(:magenta, 0.05), ylabel="Q⁻¹ ", xlabel="T(s)")

for i in eachindex(T)
    color_ = RGBf(i / 10, 0.0, 1 - i / 10)
    lines!(ax2, 1 ./ f, resp.Qinv[i, :]; color=color_, label="$(T[i] - 273)")
end
xlims!(ax2, 1e-2, 1e4)
ylims!(ax2, 10.0^(-2.5), 10.0^0.5)
fig[2, 1:2] = Legend(
    fig, ax, "Temperature (ᴼC) "; orientation=:horizontal, labelsize=12, nbanks=3)

nothing # hide
```

```@raw html
</details>
```

```@example anelastic_plts
fig # hide
```

## Master curve maxwell scaling model

```@docs; canonical = false
xfit_mxw
```

For different temperatures, the distribution with oscillation period looks like

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example anelastic_plts
T = (900:50:1400) .+ 273.0
P = 0.2
dg = 3.1
σ = 10.0 * 1.0f-3
ϕ = 0.0
ρ = 3300.0
Ch2o_ol = 0.0
T_solidus = 1100 + 273.0

f = 10 .^ -collect(range(-2, 4; length=100))

m = xfit_mxw(T, P, dg, σ, ϕ, ρ, Ch2o_ol, T_solidus, f')

resp = forward(m, []);
size(resp.J1)

fig = Figure()
ax = Axis(fig[1, 1]; xscale=log10, backgroundcolor=(:magenta, 0.05),
    ylabel="Shear Modulus (GPa)", xlabel="T(s)")

for i in eachindex(T)
    color_ = RGBf(i / 10, 0.0, 1 - i / 10)
    lines!(ax, 1 ./ f, resp.M[i, :] ./ (1e9); color=color_, label="$(T[i] - 273)")
end
xlims!(ax, 1e-2, 1e4)
ylims!(ax, 0, 80)

ax2 = Axis(fig[1, 2]; xscale=log10, yscale=log10,
    backgroundcolor=(:magenta, 0.05), ylabel="Q⁻¹ ", xlabel="T(s)")

for i in eachindex(T)
    color_ = RGBf(i / 10, 0.0, 1 - i / 10)
    lines!(ax2, 1 ./ f, resp.Qinv[i, :]; color=color_, label="$(T[i] - 273)")
end
xlims!(ax2, 1e-2, 1e4)
ylims!(ax2, 10.0^(-2.5), 10.0^0.5)
fig[2, 1:2] = Legend(
    fig, ax, "Temperature (ᴼC) "; orientation=:horizontal, labelsize=12, nbanks=2)

nothing # hide
```

```@raw html
</details>
```

```@example anelastic_plts
fig # hide
```

## Analytical Andrade model

```@docs; canonical = false
andrade_analytical
```

For different temperatures, the distribution with oscillation period looks like

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example anelastic_plts
T = (900:50:1400) .+ 273.0
P = 0.2
dg = 3.1
σ = 10.0 * 1.0f-3
ϕ = 0.0
ρ = 3300.0
Ch2o_ol = 0.0
T_solidus = 1100 + 273.0

f = 10 .^ -collect(range(-2, 4; length=100))

m = andrade_analytical(T, P, dg, σ, ϕ, ρ, Ch2o_ol, T_solidus, f')

resp = forward(m, []);
size(resp.J1)

fig = Figure()
ax = Axis(fig[1, 1]; xscale=log10, backgroundcolor=(:magenta, 0.05),
    ylabel="Shear Modulus (GPa)", xlabel="T(s)")

for i in eachindex(T)
    color_ = RGBf(i / 10, 0.0, 1 - i / 10)
    lines!(ax, 1 ./ f, resp.M[i, :] ./ (1e9); color=color_, label="$(T[i] - 273)")
end
xlims!(ax, 1e-2, 1e4)
ylims!(ax, 0, 80)

ax2 = Axis(fig[1, 2]; xscale=log10, yscale=log10,
    backgroundcolor=(:magenta, 0.05), ylabel="Q⁻¹ ", xlabel="T(s)")

for i in eachindex(T)
    color_ = RGBf(i / 10, 0.0, 1 - i / 10)
    lines!(ax2, 1 ./ f, resp.Qinv[i, :]; color=color_, label="$(T[i] - 273)")
end
xlims!(ax2, 1e-2, 1e4)
ylims!(ax2, 10.0^(-2.5), 10.0^0.5)
fig[2, 1:2] = Legend(
    fig, ax, "Temperature (ᴼC) "; orientation=:horizontal, labelsize=12, nbanks=2)

nothing # hide
```

```@raw html
</details>
```

```@example anelastic_plts
fig # hide
```
