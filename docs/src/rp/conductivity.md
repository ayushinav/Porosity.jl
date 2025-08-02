# Conductivity models

```@setup cond_plts
using MT, CairoMakie
```

## Minerals

### SEO3

```@docs; canonical = false
SEO3
```

The distribution with temperature looks like (compare with fig. 1B of [Constable, 2003](https://doi.org/10.1111/j.1365-246X.2006.03041.x)):

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example cond_plts
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

T = (700:1600) .+ 273.0
m = SEO3(T)
logsig = forward(m, []).σ

lines!(ax, inv.(T) .* 1e4, 10 .^ logsig)
lines!(ax2, inv.(T) .* 1e4, 10 .^ logsig; alpha=0)
nothing # hide
```

```@raw html
</details>
```

```@example cond_plts
f # hide
```

### Wang2006

```@docs; canonical = false
Wang2006
```

The distribution with temperature looks like (compare with fig. 2a of [Wang et al., 2006](https://www.nature.com/articles/nature05256)):

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example cond_plts
f = Figure()
ax = Axis(f[1, 1]; yscale=log10,
    xlabel="10⁴/T (K⁻¹)", ylabel="σ (S/m)",
    yticks=LogTicks(WilkinsonTicks(9; k_min=8)),
    backgroundcolor=(:magenta, 0.05))

xts = inv.([300, 600, 900, 1200, 1500] .+ 273.0) .* 1e4

ax2 = Axis(f[1, 1]; yscale=log10, xaxisposition=:top, yaxisposition=:right,
    xlabel="T (ᴼC)", xgridvisible=false, xtickformat=x -> string.(round.((1e4 ./ x) .- 273)),
    xticklabelsize=10, backgroundcolor=(:magenta, 0.05))
ax2.xticks = xts
hidespines!(ax2)
hideydecorations!(ax2)
linkyaxes!(ax, ax2)

T = (500:1400) .+ 273.0
Ch2o = [0.0, 0.01, 0.03, 0.1]' .* 1e4
m = Wang2006(T, Ch2o)
logsig = forward(m, []).σ

for i in eachindex(Ch2o)
    w = Ch2o[i]
    lines!(ax, inv.(T) .* 1e4, 10 .^ logsig[:, i]; label="$w")
    lines!(ax2, inv.(T) .* 1e4, 10 .^ logsig[:, i]; alpha=0)
end
ylims!(ax, 1e-8, 1)
ylims!(ax2, 1e-8, 1)
f[1, 2] = Legend(f, ax, "water conc. (ppm)")

nothing # hide
```

```@raw html
</details>
```

```@example cond_plts
f # hide
```

### Yoshino2009

```@docs; canonical = false
Yoshino2009
```

The distribution with temperature looks like (compare with fig. 6 of [Yoshino et al., 2009](https://doi.org/10.1016/j.epsl.2009.09.032)):

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example cond_plts
f = Figure()
ax = Axis(f[1, 1]; yscale=log10,
    xlabel="10⁴/T (K⁻¹)", ylabel="σ (S/m)",
    yticks=LogTicks(WilkinsonTicks(6; k_min=5)),
    backgroundcolor=(:magenta, 0.05))

xts = inv.([600, 800, 1000, 1200, 1400, 1600] .+ 273.0) .* 1e4

ax2 = Axis(f[1, 1]; yscale=log10, xaxisposition=:top, yaxisposition=:right,
    xlabel="T (ᴼC)", xgridvisible=false, xtickformat=x -> string.(round.((1e4 ./ x) .- 273)),
    xticklabelsize=8, backgroundcolor=(:magenta, 0.05))
ax2.xticks = xts
hidespines!(ax2)
hideydecorations!(ax2)
linkyaxes!(ax, ax2)

T = (600:1600) .+ 273.0
Ch2o = [0.0, 400, 600, 2000]'
m = Yoshino2009(T, Ch2o)
logsig = forward(m, []).σ

for i in eachindex(Ch2o)
    w = Ch2o[i]
    lines!(ax, inv.(T) .* 1e4, 10 .^ logsig[:, i]; label="$w")
    lines!(ax2, inv.(T) .* 1e4, 10 .^ logsig[:, i]; alpha=0)
end
ylims!(ax, 1e-7, 1)
ylims!(ax2, 1e-7, 1)
f[1, 2] = Legend(f, ax, "water conc. (ppm)")

nothing # hide
```

```@raw html
</details>
```

```@example cond_plts
f # hide
```

### Poe2010

```@docs; canonical = false
Poe2010
```

The distribution with temperature looks like (compare with fig. 3 of [Poe et al., 2010](https://doi.org/10.1016/j.pepi.2010.05.003)):

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example cond_plts
f = Figure()
ax = Axis(f[1, 1]; yscale=log10,
    xlabel="10⁴/T (K⁻¹)", ylabel="σ (S/m)",
    yticks=LogTicks(WilkinsonTicks(6; k_min=5)),
    backgroundcolor=(:magenta, 0.05))

xts = inv.([300, 600, 900, 1200, 1500] .+ 273.0) .* 1e4

ax2 = Axis(f[1, 1]; yscale=log10, xaxisposition=:top, yaxisposition=:right,
    xlabel="T (ᴼC)", xgridvisible=false, xtickformat=x -> string.(round.((1e4 ./ x) .- 273)),
    xticklabelsize=8, backgroundcolor=(:magenta, 0.05))
ax2.xticks = xts
hidespines!(ax2)
hideydecorations!(ax2)
linkyaxes!(ax, ax2)

T = (200:1600) .+ 273.0
Ch2o = [0.0, 400, 600, 2000]'
m = Poe2010(T, Ch2o)
logsig = forward(m, []).σ

for i in eachindex(Ch2o)
    w = Ch2o[i]
    lines!(ax, inv.(T) .* 1e4, 10 .^ logsig[:, i]; label="$w")
    lines!(ax2, inv.(T) .* 1e4, 10 .^ logsig[:, i]; alpha=0)
end
ylims!(ax, 1e-6, 10)
ylims!(ax2, 1e-6, 10)
f[1, 2] = Legend(f, ax, "water conc. (ppm)")

nothing # hide
```

```@raw html
</details>
```

```@example cond_plts
f # hide
```

### Jones2012

```@docs; canonical = false
Jones2012
```

The distribution with temperature looks like:

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example cond_plts
f = Figure()
ax = Axis(f[1, 1]; yscale=log10,
    xlabel="10⁴/T (K⁻¹)", ylabel="σ (S/m)",
    yticks=LogTicks(WilkinsonTicks(6; k_min=5)),
    backgroundcolor=(:magenta, 0.05))

xts = inv.([300, 600, 900, 1200, 1500] .+ 273.0) .* 1e4

ax2 = Axis(f[1, 1]; yscale=log10, xaxisposition=:top, yaxisposition=:right,
    xlabel="T (ᴼC)", xgridvisible=false, xtickformat=x -> string.(round.((1e4 ./ x) .- 273)),
    xticklabelsize=8, backgroundcolor=(:magenta, 0.05))
ax2.xticks = xts
hidespines!(ax2)
hideydecorations!(ax2)
linkyaxes!(ax, ax2)

T = (200:1600) .+ 273.0
Ch2o = [0.0, 400, 600, 2000]'
m = Jones2012(T, Ch2o)
logsig = forward(m, []).σ

for i in eachindex(Ch2o)
    w = Ch2o[i]
    lines!(ax, inv.(T) .* 1e4, 10 .^ logsig[:, i]; label="$w")
    lines!(ax2, inv.(T) .* 1e4, 10 .^ logsig[:, i]; alpha=0)
end
ylims!(ax, 1e-6, 10)
ylims!(ax2, 1e-6, 10)
f[1, 2] = Legend(f, ax, "water conc. (ppm)")

nothing # hide
```

```@raw html
</details>
```

```@example cond_plts
f # hide
```

### UHO2014

```@docs; canonical = false
UHO2014
```

The distribution with temperature looks like (compare with fig. 4 of [Gardés et al., 2010](https://doi.org/10.1002/2014GC005496)):

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example cond_plts
f = Figure()
ax = Axis(f[1, 1]; yscale=log10,
    xlabel="10⁴/T (K⁻¹)", ylabel="σ (S/m)",
    yticks=LogTicks(WilkinsonTicks(6; k_min=5)),
    backgroundcolor=(:magenta, 0.05))

xts = inv.([600, 800, 1000, 1200, 1400, 1600] .+ 273.0) .* 1e4

ax2 = Axis(f[1, 1]; yscale=log10, xaxisposition=:top, yaxisposition=:right,
    xlabel="T (ᴼC)", xgridvisible=false, xtickformat=x -> string.(round.((1e4 ./ x) .- 273)),
    xticklabelsize=8, backgroundcolor=(:magenta, 0.05))
ax2.xticks = xts
hidespines!(ax2)
hideydecorations!(ax2)
linkyaxes!(ax, ax2)

T = (600:1600) .+ 273.0
Ch2o = [0.0, 400, 600, 2000]'
m = Jones2012(T, Ch2o)
logsig = forward(m, []).σ

for i in eachindex(Ch2o)
    w = Ch2o[i]
    lines!(ax, inv.(T) .* 1e4, 10 .^ logsig[:, i]; label="$w")
    lines!(ax2, inv.(T) .* 1e4, 10 .^ logsig[:, i]; alpha=0)
end
ylims!(ax, 1e-10, 10)
ylims!(ax2, 1e-10, 10)
f[1, 2] = Legend(f, ax, "water conc. (ppm)")

nothing # hide
```

```@raw html
</details>
```

```@example cond_plts
f # hide
```

## Melt

  - [Sifre2014](@ref Sifre2014) : Sifre et al., 2014 : [Electrical conductivity during incipient melting in the oceanic low-velocity zone](https://doi.org/10.1038/nature13245)
  - [Gaillard2008](@ref Gaillard2008) : Gaillard et al., 2008 : [Carbonatite Melts and Electrical Conductivity in the Asthenosphere](https://www.science.org/doi/10.1126/science.1164446)

### Ni2011

```@docs; canonical = false
Ni2011
```

The distribution with temperature looks like :

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example cond_plts
f = Figure()
ax = Axis(f[1, 1]; yscale=log10,
    xlabel="10⁴/T (K⁻¹)", ylabel="σ (S/m)",
    yticks=LogTicks(WilkinsonTicks(9; k_min=8)),
    backgroundcolor=(:magenta, 0.05))

xts = inv.([1000, 1100, 1200, 1300, 1400] .+ 273.0) .* 1e4

ax2 = Axis(f[1, 1]; yscale=log10, xaxisposition=:top, yaxisposition=:right,
    xlabel="T (ᴼC)", xgridvisible=false, xtickformat=x -> string.(round.((1e4 ./ x) .- 273)),
    xticklabelsize=10, backgroundcolor=(:magenta, 0.05))
ax2.xticks = xts
hidespines!(ax2)
hideydecorations!(ax2)
linkyaxes!(ax, ax2)

T = (1000:1400) .+ 273.0
Ch2o = [0.0, 0.01, 0.03, 0.1]' .* 1e4
m = Ni2011(T, Ch2o)
logsig = forward(m, []).σ

for i in eachindex(Ch2o)
    w = Ch2o[i]
    lines!(ax, inv.(T) .* 1e4, 10 .^ logsig[:, i]; label="$w")
    lines!(ax2, inv.(T) .* 1e4, 10 .^ logsig[:, i]; alpha=0)
end

# ylims!(ax2, 1e-8, 1)
f[1, 2] = Legend(f, ax, "water conc. (ppm)")

nothing # hide
```

```@raw html
</details>
```

```@example cond_plts
f # hide
```

### Sifre2014

```@docs; canonical = false
Ni2011
```

The distribution with temperature looks like :

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example cond_plts
f = Figure()
ax = Axis(f[1, 1]; yscale=log10,
    xlabel="10⁴/T (K⁻¹)", ylabel="σ (S/m)",
    yticks=LogTicks(WilkinsonTicks(9; k_min=8)),
    backgroundcolor=(:magenta, 0.05))

xts = inv.([800, 1000, 1100, 1200, 1300, 1400] .+ 273.0) .* 1e4

ax2 = Axis(f[1, 1]; yscale=log10, xaxisposition=:top, yaxisposition=:right,
    xlabel="T (ᴼC)", xgridvisible=false, xtickformat=x -> string.(round.((1e4 ./ x) .- 273)),
    xticklabelsize=10, backgroundcolor=(:magenta, 0.05))
ax2.xticks = xts
hidespines!(ax2)
hideydecorations!(ax2)
linkyaxes!(ax, ax2)

T = (800:1400) .+ 273.0
Ch2o = [0.0, 0.01, 0.03, 0.1]' .* 1e4
Cco2_m = 1e3
m = Sifre2014(T, Ch2o, Cco2_m)
logsig = forward(m, []).σ

for i in eachindex(Ch2o)
    w = Ch2o[i]
    lines!(ax, inv.(T) .* 1e4, 10 .^ logsig[:, i]; label="$w")
    lines!(ax2, inv.(T) .* 1e4, 10 .^ logsig[:, i]; alpha=0)
end

f[1, 2] = Legend(f, ax, "water conc. (ppm)")

ax = Axis(f[2, 1]; yscale=log10,
    xlabel="10⁴/T (K⁻¹)", ylabel="σ (S/m)",
    yticks=LogTicks(WilkinsonTicks(9; k_min=8)),
    backgroundcolor=(:magenta, 0.05))

xts = inv.([800, 1000, 1100, 1200, 1300, 1400] .+ 273.0) .* 1e4

ax2 = Axis(f[2, 1]; yscale=log10, xaxisposition=:top, yaxisposition=:right,
    xlabel="T (ᴼC)", xgridvisible=false, xtickformat=x -> string.(round.((1e4 ./ x) .- 273)),
    xticklabelsize=10, backgroundcolor=(:magenta, 0.05))
ax2.xticks = xts
hidespines!(ax2)
hideydecorations!(ax2)
linkyaxes!(ax, ax2)

T = (1000:1400) .+ 273.0
Ch2o = 1e3
Cco2_m = [0.0, 0.01, 0.03, 0.1]' .* 1e4
m = Sifre2014(T, Ch2o, Cco2_m)
logsig = forward(m, []).σ

for i in eachindex(Cco2_m)
    w = Cco2_m[i]
    lines!(ax, inv.(T) .* 1e4, 10 .^ logsig[:, i]; label="$w")
    lines!(ax2, inv.(T) .* 1e4, 10 .^ logsig[:, i]; alpha=0)
end

f[2, 2] = Legend(f, ax, "CO₂ conc. (ppm)")

nothing # hide
```

```@raw html
</details>
```

```@example cond_plts
f # hide
```

### Gaillard2008

```@docs; canonical = false
Gaillard2008
```

The distribution with temperature looks like :

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example cond_plts
f = Figure()
ax = Axis(f[3:8, 3:8]; yscale=log10,
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

T = (700:1600) .+ 273.0
m = Gaillard2008(T)
logsig = forward(m, []).σ

lines!(ax, inv.(T) .* 1e4, 10 .^ logsig)
lines!(ax2, inv.(T) .* 1e4, 10 .^ logsig; alpha=0)
nothing # hide
```

```@raw html
</details>
```

```@example cond_plts
f # hide
```
