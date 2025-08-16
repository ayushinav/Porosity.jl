abstract type AbstractViscousModel <: AbstractRockphyModel end

mutable struct RockphyViscous{T1, T2} <: AbstractRockphyResponse
    ϵ_rate::T1
    η::T2
end

"""
    HZK2011(T, P, dg, σ, ϕ)

Calculate strain rate and viscosity for steady state olivine flow,
per Zimmerman and Kohlstedt (2011), using the three creep mechanisms, i.e.,
diffusion, dislocation, grain boundary sliding

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `dg`: Grain size (μm)
    - `σ` : Shear stress (GPa)
    - `ϕ` : Porosity

## Keyword Arguments

    - `params` : Various coefficients required for calculation.
    Coefficients for different mechanisms (stored in `mechs` field):
        - `diff` : Diffusion creep
        - `disl` : Dislocation creep
        - `gbs`  : Grain boundary sliding 
    To investigate coefficients, call `default_params(Val{HZK2011})`. 
    To modify coefficients, check the relevant documentation page. This
    will also users to get any particular type of mechanism, eg. `diff` only
    by setting the `A` in `disl` and `gbs` to `0f0`.

    `params` for `HZK2011` holds another important field:
        - `melt_enhancement` : TODO

## Usage

```julia
T = collect(1073.0f0:30:1373.0f0)
P = 2 .+ zero(T)
dg = collect(3.0f0:4.0f-1:7.0f0)
σ = collect(7.5f0:0.5f0:12.5f0) .* 1.0f-3
ϕ = collect(1.0f-2:1.0f-3:2.0f-2)

model = HZK2011(T, P, dg, σ, ϕ)

forward(model, [])
```

## References

  - Hansen, Zimmerman and Kohlstedt, 2011, "Grain boundary sliding in San Carlos olivine:
    Flow law parameters and crystallographic-preferred orientation", J. Geophys. Res.,
    https://doi.org/10.1029/2011JB008220
"""
mutable struct HZK2011{T1, T2, T3, T4, T5} <: AbstractViscousModel
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
end

"""
    HK2003(T, P, dg, σ, ϕ, Ch2o_ol = 0.)

Calculate strain rate and viscosity for steady state olivine flow,
per Hirth and Kohlstedt (2003), using the three creep mechanisms, i.e.,
diffusion, dislocation, grain boundary sliding

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `dg`: Grain size (μm)
    - `σ` : Shear stress (GPa)
    - `ϕ` : Porosity

## Optional Arguments

    - `Ch2o_ol` : Water concentration in olivine (ppm), defaults to 0 ppm.

## Keyword Arguments

    - `params` : Various coefficients required for calculation.
    Coefficients for different mechanisms (stored in `mechs` field):
        - `diff` : Diffusion creep
        - `disl` : Dislocation creep
        - `gbs`  : Grain boundary sliding 
    To investigate coefficients, call `default_params(Val{HK2003})`. 
    To modify coefficients, check the relevant documentation page. This
    will also users to get any particular type of mechanism, eg. `diff` only
    by setting the `A` in `disl` and `gbs` to `0f0`.

    `params` for `HK2003` holds another important field:
    - `melt_enhancement` : TODO

## Usage

```julia
T = collect(1073.0f0:30:1373.0f0)
P = 2 .+ zero(T)
dg = collect(3.0f0:4.0f-1:7.0f0)
σ = collect(7.5f0:0.5f0:12.5f0) .* 1.0f-3
ϕ = collect(1.0f-2:1.0f-3:2.0f-2)

model = HK2003(T, P, dg, σ, ϕ)

forward(model, [])
```

## References

  - Hirth and Kohlstedt, 2003, "Rheology of the Upper Mantle and the Mantle Wedge: A View from the Experimentalists",
    Inside the Subduction Factory, J. Eiler (Ed.).
    https://doi.org/10.1029/138GM06
"""
mutable struct HK2003{T1, T2, T3, T4, T5, T6} <: AbstractViscousModel
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    Ch2o_ol::T6
end

HK2003(T, P, dg, σ, ϕ) = HK2003(T, P, dg, σ, ϕ, 0.0f0)

"""
    xfit_premelt(T, P, dg, σ, ϕ)

Calculate viscosity for steady state olivine flow for pre-melting, i.e.,
temperatures are just below and above the solidus, per Yamauchi and Takei (2016)

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `dg`: Grain size (μm)
    - `σ` : Shear stress (GPa)
    - `ϕ` : Porosity
    - `T_solidus` : Solidus temperature

## Keyword Arguments

    - `params` : Various coefficients required for calculation.
    Coefficients for different mechanisms (stored in `mechs` field):
        - `diff` : Diffusion creep
        - `disl` : Dislocation creep
        - `gbs`  : Grain boundary sliding 
    To investigate coefficients, call `default_params(Val{xfit_premelt}())`. 
    To modify coefficients, check the relevant documentation page. This
    will also users to get any particular type of mechanism, eg. `diff` only
    by setting the `A` in `disl` and `gbs` to `0f0`.

## Usage

```julia
T = collect(1073.0f0:30:1373.0f0)
P = 2 .+ zero(T)
dg = collect(3.0f0:4.0f-1:7.0f0)
σ = collect(7.5f0:0.5f0:12.5f0) .* 1.0f-3
ϕ = collect(1.0f-2:1.0f-3:2.0f-2)
T_solidus = 1200 + 273 .+ zero(T)

model = xfit_premelt(T, P, dg, σ, ϕ, T_solidus)

forward(model, [])
```

## References

  - Yamauchi and Takei, 2016, "Polycrystal anelasticity at near-solidus temperatures",
    J. Geophys. Res. Solid Earth,
    https://doi.org/10.1002/2016JB013316
"""
mutable struct xfit_premelt{T1, T2, T3, T4, T5, T6} <: AbstractViscousModel
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    T_solidus::T6
end

default_params(::Type{T}) where {T <: HZK2011} = default_params_HZK2011
default_params(::Type{T}) where {T <: HK2003} = default_params_HK2003
default_params(::Type{T}) where {T <: xfit_premelt} = default_params_xfit_premelt
