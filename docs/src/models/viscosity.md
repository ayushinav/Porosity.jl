# Viscosity models

```@setup visc_plts
using Porosity, CairoMakie
```

## HK2003

```@docs; canonical = false
HK2003
```

### Varying grain size

The distribution with temperature looks like (assuming increase in grain size is independent of other parameters) at 2 GPa pressure with a porosity of 0.015:

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example visc_plts
f = Figure(; size=(700, 400))

T = (600:1600) .+ 273.0
P = 2.0
dg = [1.0, 3.0, 10.0, 30.0]'
σ = 10e-3
ϕ = 0.015

m = HK2003(T, P, dg, σ, ϕ)
resp = forward(m, []);

resp_fields = [:ϵ_rate, :η]
units = ["", "(Pa s)"]

ax_coords = [(1, 1), (1, 2)]
for i in eachindex(resp_fields)
    ax = Axis(f[ax_coords[i]...];
        xlabel="10⁴/T (K⁻¹)", ylabel=string(resp_fields[i]) * " $(units[i])",
        yticks=LogTicks(WilkinsonTicks(6; k_min=5)),
        backgroundcolor=(:magenta, 0.052))

    xts = inv.([600, 800, 1000, 1200, 1600] .+ 273.0) .* 1e4

    ax2 = Axis(f[ax_coords[i]...]; xaxisposition=:top, yaxisposition=:right,
        xlabel="T (ᴼC)", xgridvisible=false,
        xtickformat=x -> string.(round.((1e4 ./ x) .- 273)), xticklabelsize=8,
        backgroundcolor=(:magenta, 0.05))
    ax2.xticks = xts
    hidespines!(ax2)
    hideydecorations!(ax2)
    linkyaxes!(ax, ax2)
    ylims!(ax, extrema(getfield(resp, resp_fields[i])) .* (0.3, 3))
    ylims!(ax2, extrema(getfield(resp, resp_fields[i])) .* (0.3, 3))
    ax.yscale = log10
    ax2.yscale = log10

    for j in axes(getfield(resp, resp_fields[i]), 2)
        d = dg[j]
        r_ = getfield(resp, resp_fields[i])
        lines!(ax, inv.(T) .* 1e4, r_[:, j]; label="$d")
        lines!(ax2, inv.(T) .* 1e4, r_[:, j]; alpha=0)
    end
end
f[2, 1:2] = Legend(f, f.content[end - 1], "grain size (μm)"; orientation=:horizontal)

nothing # hide
```

```@raw html
</details>
```

```@example visc_plts
f # hide
```

### Varying porosity

We can also look at the distribution with parameters by varying just the porosity :

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example visc_plts
f = Figure(; size=(700, 400))

T = (600:1600) .+ 273.0
P = 2.0
dg = 4.0
σ = 10e-3
ϕ = [0.01, 0.03, 0.1, 0.3]'

m = HK2003(T, P, dg, σ, ϕ)
resp = forward(m, []);

resp_fields = [:ϵ_rate, :η]
units = ["", "(Pa⋅s)"]

ax_coords = [(1, 1), (1, 2)]
for i in eachindex(resp_fields)
    ax = Axis(f[ax_coords[i]...];
        xlabel="10⁴/T (K⁻¹)", ylabel=string(resp_fields[i]) * " $(units[i])",
        yticks=LogTicks(WilkinsonTicks(6; k_min=5)),
        backgroundcolor=(:magenta, 0.052))

    xts = inv.([600, 800, 1000, 1200, 1600] .+ 273.0) .* 1e4

    ax2 = Axis(f[ax_coords[i]...]; xaxisposition=:top, yaxisposition=:right,
        xlabel="T (ᴼC)", xgridvisible=false,
        xtickformat=x -> string.(round.((1e4 ./ x) .- 273)), xticklabelsize=8,
        backgroundcolor=(:magenta, 0.05))
    ax2.xticks = xts
    hidespines!(ax2)
    hideydecorations!(ax2)
    linkyaxes!(ax, ax2)
    ylims!(ax, extrema(getfield(resp, resp_fields[i])) .* (0.3, 3))
    ylims!(ax2, extrema(getfield(resp, resp_fields[i])) .* (0.3, 3))
    ax.yscale = log10
    ax2.yscale = log10

    for j in axes(getfield(resp, resp_fields[i]), 2)
        d = ϕ[j]
        r_ = getfield(resp, resp_fields[i])
        lines!(ax, inv.(T) .* 1e4, r_[:, j]; label="$d")
        lines!(ax2, inv.(T) .* 1e4, r_[:, j]; alpha=0)
    end
end
f[2, 1:2] = Legend(f, f.content[end - 1], "porosity"; orientation=:horizontal)
nothing # hide
```

```@raw html
</details>
```

```@example visc_plts
f # hide
```

### Varying water content

Another important parameter to consider for `HK2003` is water content. By default, we assume 0 ppm. The distribution with temperature for different water content (keeping other parameters constant) looks like :

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example visc_plts
f = Figure(; size=(700, 400))

T = (600:1600) .+ 273.0
P = 2.0
dg = 4.0
σ = 10e-3
ϕ = 0.015
Ch2o = [0.0, 100.0, 300.0, 1000.0]'

m = HK2003(T, P, dg, σ, ϕ, Ch2o)
resp = forward(m, []);

resp_fields = [:ϵ_rate, :η]
units = ["", "(Pa⋅s)"]

ax_coords = [(1, 1), (1, 2)]
for i in eachindex(resp_fields)
    ax = Axis(f[ax_coords[i]...];
        xlabel="10⁴/T (K⁻¹)", ylabel=string(resp_fields[i]) * " $(units[i])",
        yticks=LogTicks(WilkinsonTicks(6; k_min=5)),
        backgroundcolor=(:magenta, 0.052))

    xts = inv.([600, 800, 1000, 1200, 1600] .+ 273.0) .* 1e4

    ax2 = Axis(f[ax_coords[i]...]; xaxisposition=:top, yaxisposition=:right,
        xlabel="T (ᴼC)", xgridvisible=false,
        xtickformat=x -> string.(round.((1e4 ./ x) .- 273)), xticklabelsize=8,
        backgroundcolor=(:magenta, 0.05))
    ax2.xticks = xts
    hidespines!(ax2)
    hideydecorations!(ax2)
    linkyaxes!(ax, ax2)
    ylims!(ax, extrema(getfield(resp, resp_fields[i])) .* (0.3, 3))
    ylims!(ax2, extrema(getfield(resp, resp_fields[i])) .* (0.3, 3))
    ax.yscale = log10
    ax2.yscale = log10

    for j in axes(getfield(resp, resp_fields[i]), 2)
        d = Ch2o[j]
        r_ = getfield(resp, resp_fields[i])
        lines!(ax, inv.(T) .* 1e4, r_[:, j]; label="$d")
        lines!(ax2, inv.(T) .* 1e4, r_[:, j]; alpha=0)
    end
end
f[2, 1:2] = Legend(f, f.content[end - 1], "Water conc. (ppm)"; orientation=:horizontal)
nothing # hide
```

```@raw html
</details>
```

```@example visc_plts
f # hide
```

## HZK2011

```@docs; canonical = false
HZK2011
```

### Varying grain size

The distribution with temperature looks like (assuming increase in grain size is independent of other parameters) at 2 GPa pressure with a porosity of 0.015:

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example visc_plts
f = Figure(; size=(700, 400))

T = (600:1600) .+ 273.0
P = 2.0
dg = [1.0, 3.0, 10.0, 30.0]'
σ = 10e-3
ϕ = 0.015

m = HZK2011(T, P, dg, σ, ϕ)
resp = forward(m, []);

resp_fields = [:ϵ_rate, :η]
units = ["", "(Pa s)"]

ax_coords = [(1, 1), (1, 2)]
for i in eachindex(resp_fields)
    ax = Axis(f[ax_coords[i]...];
        xlabel="10⁴/T (K⁻¹)", ylabel=string(resp_fields[i]) * " $(units[i])",
        yticks=LogTicks(WilkinsonTicks(6; k_min=5)),
        backgroundcolor=(:magenta, 0.052))

    xts = inv.([600, 800, 1000, 1200, 1600] .+ 273.0) .* 1e4

    ax2 = Axis(f[ax_coords[i]...]; xaxisposition=:top, yaxisposition=:right,
        xlabel="T (ᴼC)", xgridvisible=false,
        xtickformat=x -> string.(round.((1e4 ./ x) .- 273)), xticklabelsize=8,
        backgroundcolor=(:magenta, 0.05))
    ax2.xticks = xts
    hidespines!(ax2)
    hideydecorations!(ax2)
    linkyaxes!(ax, ax2)
    ylims!(ax, extrema(getfield(resp, resp_fields[i])) .* (0.3, 3))
    ylims!(ax2, extrema(getfield(resp, resp_fields[i])) .* (0.3, 3))
    ax.yscale = log10
    ax2.yscale = log10

    for j in axes(getfield(resp, resp_fields[i]), 2)
        d = dg[j]
        r_ = getfield(resp, resp_fields[i])
        lines!(ax, inv.(T) .* 1e4, r_[:, j]; label="$d")
        lines!(ax2, inv.(T) .* 1e4, r_[:, j]; alpha=0)
    end
end
f[2, 1:2] = Legend(f, f.content[end - 1], "grain size (μm)"; orientation=:horizontal)

nothing # hide
```

```@raw html
</details>
```

```@example visc_plts
f # hide
```

### Varying porosity

We can also look at the distribution with parameters by varying just the porosity :

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example visc_plts
f = Figure(; size=(700, 400))

T = (600:1600) .+ 273.0
P = 2.0
dg = 4.0
σ = 10e-3
ϕ = [0.01, 0.03, 0.1, 0.3]'

m = HZK2011(T, P, dg, σ, ϕ)
resp = forward(m, []);

resp_fields = [:ϵ_rate, :η]
units = ["", "(Pa⋅s)"]

ax_coords = [(1, 1), (1, 2)]
for i in eachindex(resp_fields)
    ax = Axis(f[ax_coords[i]...];
        xlabel="10⁴/T (K⁻¹)", ylabel=string(resp_fields[i]) * " $(units[i])",
        yticks=LogTicks(WilkinsonTicks(6; k_min=5)),
        backgroundcolor=(:magenta, 0.052))

    xts = inv.([600, 800, 1000, 1200, 1600] .+ 273.0) .* 1e4

    ax2 = Axis(f[ax_coords[i]...]; xaxisposition=:top, yaxisposition=:right,
        xlabel="T (ᴼC)", xgridvisible=false,
        xtickformat=x -> string.(round.((1e4 ./ x) .- 273)), xticklabelsize=8,
        backgroundcolor=(:magenta, 0.05))
    ax2.xticks = xts
    hidespines!(ax2)
    hideydecorations!(ax2)
    linkyaxes!(ax, ax2)
    ylims!(ax, extrema(getfield(resp, resp_fields[i])) .* (0.3, 3))
    ylims!(ax2, extrema(getfield(resp, resp_fields[i])) .* (0.3, 3))
    ax.yscale = log10
    ax2.yscale = log10

    for j in axes(getfield(resp, resp_fields[i]), 2)
        d = ϕ[j]
        r_ = getfield(resp, resp_fields[i])
        lines!(ax, inv.(T) .* 1e4, r_[:, j]; label="$d")
        lines!(ax2, inv.(T) .* 1e4, r_[:, j]; alpha=0)
    end
end
f[2, 1:2] = Legend(f, f.content[end - 1], "porosity"; orientation=:horizontal)
nothing # hide
```

```@raw html
</details>
```

```@example visc_plts
f # hide
```

## xfit_premelt

```@docs; canonical = false
xfit_premelt
```

!!!note
Note that `xfit_premelt` populates the return value of strain rate with zeros. We, therefore, do not plot the same here.

### Varying grain size

The distribution with temperature looks like (assuming increase in grain size is independent of other parameters) at 2 GPa pressure with a porosity of 0.015 at solidus temperature of 1200 ᴼC:

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example visc_plts
f = Figure(; size=(600, 400))

T = (600:1600) .+ 273.0
P = 2.0
dg = [1.0, 3.0, 10.0, 30.0]'
σ = 10e-3
ϕ = 0.015
T_solidus = 1200 + 273

m = xfit_premelt(T, P, dg, σ, ϕ, T_solidus)
resp = forward(m, []);

resp_fields = [:η]
units = ["(Pa s)"]

ax_coords = [(1, 1), (1, 2)]
for i in eachindex(resp_fields)
    ax = Axis(f[ax_coords[i]...];
        xlabel="10⁴/T (K⁻¹)", ylabel=string(resp_fields[i]) * " $(units[i])",
        yticks=LogTicks(WilkinsonTicks(6; k_min=5)),
        backgroundcolor=(:magenta, 0.052))

    xts = inv.([600, 800, 1000, 1200, 1600] .+ 273.0) .* 1e4

    ax2 = Axis(f[ax_coords[i]...]; xaxisposition=:top, yaxisposition=:right,
        xlabel="T (ᴼC)", xgridvisible=false,
        xtickformat=x -> string.(round.((1e4 ./ x) .- 273)), xticklabelsize=8,
        backgroundcolor=(:magenta, 0.05))
    ax2.xticks = xts
    hidespines!(ax2)
    hideydecorations!(ax2)
    linkyaxes!(ax, ax2)

    if sum(extrema(getfield(resp, resp_fields[i]))) > 1e-15
        ylims!(ax, extrema(getfield(resp, resp_fields[i])) .* (0.3, 3))
        ylims!(ax2, extrema(getfield(resp, resp_fields[i])) .* (0.3, 3))
        ax.yscale = log10
        ax2.yscale = log10
    end
    for j in axes(getfield(resp, resp_fields[i]), 2)
        d = dg[j]
        r_ = getfield(resp, resp_fields[i])
        lines!(ax, inv.(T) .* 1e4, r_[:, j]; label="$d")
        lines!(ax2, inv.(T) .* 1e4, r_[:, j]; alpha=0)
    end
end

f[1, 2] = Legend(f, f.content[end - 1], "grain size (μm)")

for i in eachindex(f.content[1:(end - 1)])
    ax_ = f.content[i]
    lv = vlines!(ax_, 1e4 / T_solidus; color=:green, linestyle=:dash)
    axislegend(
        ax_, [lv], ["solidus temperature"]; position=:rb, labelsize=12, framevisible=false)
end

nothing # hide
```

```@raw html
</details>
```

```@example visc_plts
f # hide
```

### Varying porosity

We can also look at the distribution with parameters by varying just the porosity :

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example visc_plts
f = Figure(; size=(600, 400))

T = (600:1600) .+ 273.0
P = 2.0
dg = 4.0
σ = 10e-3
ϕ = [0.01, 0.03, 0.1, 0.3]'
T_solidus = 1200 + 273

m = xfit_premelt(T, P, dg, σ, ϕ, T_solidus)
resp = forward(m, []);

resp_fields = [:η]
units = ["(Pa s)"]

ax_coords = [(1, 1), (1, 2)]
for i in eachindex(resp_fields)
    ax = Axis(f[ax_coords[i]...];
        xlabel="10⁴/T (K⁻¹)", ylabel=string(resp_fields[i]) * " $(units[i])",
        yticks=LogTicks(WilkinsonTicks(6; k_min=5)),
        backgroundcolor=(:magenta, 0.052))

    xts = inv.([600, 800, 1000, 1200, 1600] .+ 273.0) .* 1e4

    ax2 = Axis(f[ax_coords[i]...]; xaxisposition=:top, yaxisposition=:right,
        xlabel="T (ᴼC)", xgridvisible=false,
        xtickformat=x -> string.(round.((1e4 ./ x) .- 273)), xticklabelsize=8,
        backgroundcolor=(:magenta, 0.05))
    ax2.xticks = xts
    hidespines!(ax2)
    hideydecorations!(ax2)
    linkyaxes!(ax, ax2)

    if sum(extrema(getfield(resp, resp_fields[i]))) > 1e-15
        ylims!(ax, extrema(getfield(resp, resp_fields[i])) .* (0.3, 3))
        ylims!(ax2, extrema(getfield(resp, resp_fields[i])) .* (0.3, 3))
        ax.yscale = log10
        ax2.yscale = log10
    end
    for j in axes(getfield(resp, resp_fields[i]), 2)
        d = ϕ[j]
        r_ = getfield(resp, resp_fields[i])
        lines!(ax, inv.(T) .* 1e4, r_[:, j]; label="$d")
        lines!(ax2, inv.(T) .* 1e4, r_[:, j]; alpha=0)
    end
end

f[1, 2] = Legend(f, f.content[end - 1], "Porosity")

for i in eachindex(f.content[1:(end - 1)])
    ax_ = f.content[i]
    lv = vlines!(ax_, 1e4 / T_solidus; color=:green, linestyle=:dash)
    axislegend(
        ax_, [lv], ["solidus temperature"]; position=:rb, labelsize=12, framevisible=false)
end

nothing # hide
```

```@raw html
</details>
```

```@example visc_plts
f # hide
```
