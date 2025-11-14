# Elasticity models

```@setup elastic_plts
using Porosity, CairoMakie
```

## Anharmonic scaling

```@docs; canonical = false
anharmonic
```

The distribution with temperature looks like (assuming increase in density is not due to increase in pressure) at 2 GPa pressure:

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example elastic_plts
f = Figure(; size=(700, 700))

T = (600:1600) .+ 273.0
P = 2.0
ρ = [2700, 2900, 3100, 3300]'
m = anharmonic(T, P, ρ)
resp = forward(m, []);

resp_fields = [:K, :G, :Vp, :Vs]
units = ["GPa", "GPa", "km/s", "km/s"]

ax_coords = [(1, 1), (1, 2), (2, 1), (2, 2)]
for i in eachindex(resp_fields)
    ax = Axis(f[ax_coords[i]...];
        xlabel="10⁴/T (K⁻¹)", ylabel=string(resp_fields[i]) * " ($(units[i]))",
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

    for j in axes(getfield(resp, resp_fields[i]), 2)
        d = ρ[j]
        r_ = getfield(resp, resp_fields[i])
        lines!(ax, inv.(T) .* 1e4, r_[:, j]; label="$d")
        lines!(ax2, inv.(T) .* 1e4, r_[:, j]; alpha=0)
    end
end
f[3, 1:2] = Legend(f, f.content[end - 1], "density (kg/m³)"; orientation=:horizontal)

nothing # hide
```

```@raw html
</details>
```

```@example elastic_plts
f # hide
```

## Anharmonic (porous) scaling

```@docs; canonical = false
anharmonic
```

The distribution with temperature looks like (assuming increase in porosity is not related to pressure (2 GPa) or density (3000 kg/m³)) :

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example elastic_plts
f = Figure(; size=(700, 700))

T = (600:1600) .+ 273.0
P = 2.0
ρ = 3300.0
ϕ = [0.01, 0.03, 0.1, 0.3]'
m = anharmonic_poro(T, P, ρ, ϕ)
resp = forward(m, []);

resp_fields = [:K, :G, :Vp, :Vs]
units = ["GPa", "GPa", "km/s", "km/s"]

ax_coords = [(1, 1), (1, 2), (2, 1), (2, 2)]
for i in eachindex(resp_fields)
    ax = Axis(f[ax_coords[i]...];
        xlabel="10⁴/T (K⁻¹)", ylabel=string(resp_fields[i]) * " ($(units[i]))",
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

    for j in axes(getfield(resp, resp_fields[i]), 2)
        d = ϕ[j]
        r_ = getfield(resp, resp_fields[i])
        lines!(ax, inv.(T) .* 1e4, r_[:, j]; label="$d")
        lines!(ax2, inv.(T) .* 1e4, r_[:, j]; alpha=0)
    end
end
f[3, 1:2] = Legend(f, f.content[end - 1], "Porosity"; orientation=:horizontal)

nothing # hide
```

```@raw html
</details>
```

```@example elastic_plts
f # hide
```

The distribution with temperature at constant porosity (0.1) (assuming increase in density is not related to pressure (2 GPa)) :

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example elastic_plts
f = Figure(; size=(700, 700))

T = (600:1600) .+ 273.0
P = 2.0
ρ = [2700, 2900, 3100, 3300]'
ϕ = 0.1
m = anharmonic_poro(T, P, ρ, ϕ)
resp = forward(m, []);

resp_fields = [:K, :G, :Vp, :Vs]
units = ["GPa", "GPa", "km/s", "km/s"]

ax_coords = [(1, 1), (1, 2), (2, 1), (2, 2)]
for i in eachindex(resp_fields)
    ax = Axis(f[ax_coords[i]...];
        xlabel="10⁴/T (K⁻¹)", ylabel=string(resp_fields[i]) * " ($(units[i]))",
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

    for j in axes(getfield(resp, resp_fields[i]), 2)
        d = ρ[j]
        r_ = getfield(resp, resp_fields[i])
        lines!(ax, inv.(T) .* 1e4, r_[:, j]; label="$d")
        lines!(ax2, inv.(T) .* 1e4, r_[:, j]; alpha=0)
    end
end
f[3, 1:2] = Legend(f, f.content[end - 1], "Porosity"; orientation=:horizontal)

nothing # hide
```

```@raw html
</details>
```

```@example elastic_plts
f # hide
```

## SLB2005

```@docs; canonical = false
SLB2005
```

The distribution with temperature looks like (assuming increase in density is not due to increase in pressure) at 2 GPa pressure:

!!! note
    
    Note that the `K`, `G` and `Vp` fields are populated with zeros.

```@raw html
<details closed><summary>Code for this figure</summary>
```

```@example elastic_plts
f = Figure(; size=(700, 700))

T = (600:1600) .+ 273.0
P = 2.0
m = SLB2005(T, P)
resp = forward(m, []);

resp_fields = [:K, :G, :Vp, :Vs]
units = ["GPa", "GPa", "km/s", "km/s"]

ax_coords = [(1, 1), (1, 2), (2, 1), (2, 2)]
for i in eachindex(resp_fields)
    ax = Axis(f[ax_coords[i]...];
        xlabel="10⁴/T (K⁻¹)", ylabel=string(resp_fields[i]) * " ($(units[i]))",
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

    r_ = getfield(resp, resp_fields[i])
    lines!(ax, inv.(T) .* 1e4, r_)
    lines!(ax2, inv.(T) .* 1e4, r_; alpha=0)
end

nothing # hide
```

```@raw html
</details>
```

```@example elastic_plts
f # hide
```
