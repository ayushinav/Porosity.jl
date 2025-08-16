# Mixing phases

```@setup mix_phases
using Porosity, CairoMakie
```

For various rock physics models, we require the bulk property by combing multiple phases. We allow the feature to mix multiple phases to get the bulk property.

!!! todo
    
    Bulk electrical conductivity of two phases is only supported at the moment.

Hashin-Shtrikman bounds, provided through [`HS1962_plus`](@ref) and [`HS1962_minus`](@ref), are often used to estimate the bulk conductivity when two phases are mixed. For conductivity, modified Archie's law is also provided through [`MAL`](@ref).

To get things started, we first need to define the phases we need to mix the models and the mixing law. This is conveniently done using [`two_phase_modelType`](@ref)

```@example mix_phases
m = two_phase_modelType(Yoshino2009, Sifre2014, HS1962_plus())
```

Now, we define the parameters required to define the model.

```@example mix_phases
T = collect(1000:1400) .+ 273.0f0
Ch2o_ol = 1.0f0
Ch2o_m = 1000.0f0
Cco2_m = 10.0f0
ϕ = 0.1f0

ps_nt = (; ϕ=ϕ, T=T, Ch2o_ol=Ch2o_ol, Ch2o_m=Ch2o_m, Cco2_m=Cco2_m)
model = m(ps_nt)
nothing # hide
```

and then as usual get the response

```@example mix_phases
resp = forward(model, [])
nothing # hide
```

and then plot :

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example mix_phases
f = Figure()
ax = Axis(f[1, 1]; yscale=log10,
    xlabel="10⁴/T (K⁻¹)", ylabel="σ (S/m)",
    yticks=LogTicks(WilkinsonTicks(6; k_min=5)),
    backgroundcolor=(:magenta, 0.05))

xts = inv.([700, 900, 1100, 1300, 1500] .+ 273.0) .* 1e4

ax2 = Axis(f[1, 1]; yscale=log10, xaxisposition=:top, yaxisposition=:right,
    xlabel="T (ᴼC)", xgridvisible=false, xtickformat=x -> string.(round.((1e4 ./ x) .- 273)),
    xticklabelsize=10, backgroundcolor=(:magenta, 0.05))
ax2.xticks = xts
hidespines!(ax2)
hideydecorations!(ax2)
linkyaxes!(ax, ax2)

logsig = resp.σ

lines!(ax, inv.(T) .* 1e4, 10 .^ logsig)
lines!(ax2, inv.(T) .* 1e4, 10 .^ logsig; alpha=0)
nothing # hide
```

```@raw html
</details>
```

```@example mix_phases
f # hide
```

The above might seem complicated at the first sight but bears a very close analogy with the other rock physics types. Lets break it down bottom up :

  - We had `forward(model, [])` similar to any other rock physics type, the last step to get the responses is to call the `forward` function.
  - Before that, we had `m(ps_nt)` where `ps_nt` had the parameters. This looks very much like `SEO3(1000. + 273)` or `anharmonic(T, P, ρ)`. Instead of passing the parameters as such, we now have to pass them through the `NamedTuple` called `ps_nt` because now, we do not know beforehand what parameters the mixed model will depend on. Had we used `SEO3` with `Gaillard2008`, we would have only needed temperature (and obviously melt fraction).
  - Before that, in the very first step, we had `two_phase_modelType(Yoshino2009, Sifre2014, HS1962_plus())`. Now, we do not have a comparison step here but we do know that the output from this `m` is used in the similar fashion as `SEO3`, `Yoshino2009` or `Sifre2014`. The `two_phase_modelType` function allows us to create one of these "types". For two phase mixing, we require the two phases along with the mixing law. Once we know them, we can completely define the physics at play, similar to how `SEO3` or `Yoshino2009` defines it in their way.

For different models, we have the distribution with melt fraction as (also a nice example of broadcasting) :

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example mix_phases
T = collect(1000:100:1400) .+ 273.0f0
Ch2o_ol = 1.0f0
Ch2o_m = 1000.0f0
Cco2_m = 10.0f0
ϕ = (0.01f0:0.01f0:0.2f0)'

m = two_phase_modelType(Yoshino2009, Sifre2014, HS1962_plus())
ps_nt = (; ϕ=ϕ, T=T, Ch2o_ol=Ch2o_ol, Ch2o_m=Ch2o_m, Cco2_m=Cco2_m)
model = m(ps_nt)

sig1 = 10.0f0 .^ forward(model, []).σ

m = two_phase_modelType(Yoshino2009, Sifre2014, HS1962_minus())
ps_nt = (; ϕ=ϕ, T=T, Ch2o_ol=Ch2o_ol, Ch2o_m=Ch2o_m, Cco2_m=Cco2_m)
model = m(ps_nt)

sig2 = 10.0f0 .^ forward(model, []).σ

m = two_phase_modelType(Yoshino2009, Sifre2014, MAL(0.2f0))
ps_nt = (; ϕ=ϕ, T=T, Ch2o_ol=Ch2o_ol, Ch2o_m=Ch2o_m, Cco2_m=Cco2_m)
model = m(ps_nt)

sig3 = 10.0f0 .^ forward(model, []).σ

# plots

fig = Figure(; size=(900, 400))
ax1 = Axis(fig[1, 1]; yscale=log10, backgroundcolor=(:magenta, 0.05),
    ylabel="σ (S/m)", xlabel="ϕ", title="HS+ bound")

for i in eachindex(T)
    color_ = RGBf(i / 5, 0.0, 1 - i / 5)
    lines!(ax1, ϕ[:], sig1[i, :]; color=color_, label="$(T[i] - 273)")
end

ax2 = Axis(fig[1, 2]; yscale=log10, backgroundcolor=(:magenta, 0.05),
    ylabel="σ (S/m)", xlabel="T(s)", title="HS- bound")

for i in eachindex(T)
    color_ = RGBf(i / 5, 0.0, 1 - i / 5)
    lines!(ax2, ϕ[:], sig2[i, :]; color=color_, label="$(T[i] - 273)")
end

ax3 = Axis(fig[1, 3]; yscale=log10, backgroundcolor=(:magenta, 0.05),
    ylabel="σ (S/m)", xlabel="T(s)", title="Modified Archie's law")

for i in eachindex(T)
    color_ = RGBf(i / 5, 0.0, 1 - i / 5)
    lines!(ax3, ϕ[:], sig3[i, :]; color=color_, label="$(T[i] - 273)")
end

fig[2, 1:3] = Legend(
    fig, ax3, "Temperature (ᴼC)"; orientation=:horizontal, labelsize=12)
nothing # hide
```

```@raw html
</details>
```

```@example mix_phases
fig # hide
```
