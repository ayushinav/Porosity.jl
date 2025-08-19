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

```jldoctest
julia> T = collect(1073.0f0:40:1273.0f0);

julia> P = 2 .+ zero(T);

julia> dg = collect(3.0f0:8.0f-1:7.0f0);

julia> σ = collect(7.5f0:1.0f0:12.5f0) .* 1.0f-3;

julia> ϕ = collect(1.0f-2:2.0f-3:2.0f-2);

julia> model = HZK2011(T, P, dg, σ, ϕ)
Model : HZK2011
Temperature (K) : Float32[1073.0, 1113.0, 1153.0, 1193.0, 1233.0, 1273.0]
Pressure (GPa) : Float32[2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
grain size(μm) : Float32[3.0, 3.8, 4.6, 5.4, 6.2, 7.0]
Shear stress (GPa) : Float32[0.0075000003, 0.0085, 0.009500001, 0.010500001, 0.011500001, 0.0125]
Porosity : Float32[0.01, 0.012, 0.014, 0.016, 0.018, 0.02]

julia> forward(model, [])
Rock physics Response : RockphyViscous
Strain rate : Float32[8.3687815f-13, 2.4096476f-12, 7.02197f-12, 2.0106366f-11, 5.5824817f-11, 1.4951564f-10]
Viscosity (Pa s) : Float32[8.961879f18, 3.5274867f18, 1.3528968f18, 5.222227f17, 2.060016f17, 8.36033f16]
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

```jldoctest
julia> T = collect(1073.0f0:40:1273.0f0);

julia> P = 2 .+ zero(T);

julia> dg = collect(3.0f0:8.0f-1:7.0f0);

julia> σ = collect(7.5f0:1.0f0:12.5f0) .* 1.0f-3;

julia> ϕ = collect(1.0f-2:2.0f-3:2.0f-2);

julia> model = HK2003(T, P, dg, σ, ϕ)
Model : HK2003
Temperature (K) : Float32[1073.0, 1113.0, 1153.0, 1193.0, 1233.0, 1273.0]
Pressure (GPa) : Float32[2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
grain size(μm) : Float32[3.0, 3.8, 4.6, 5.4, 6.2, 7.0]
Shear stress (GPa) : Float32[0.0075000003, 0.0085, 0.009500001, 0.010500001, 0.011500001, 0.0125]
Porosity : Float32[0.01, 0.012, 0.014, 0.016, 0.018, 0.02]
Water concentration in olivine (ppm) : 0.0

julia> forward(model, [])
Rock physics Response : RockphyViscous
Strain rate : Float32[3.152961f-11, 9.076501f-11, 2.643556f-10, 7.560147f-10, 2.093749f-9, 5.5805227f-9]
Viscosity (Pa s) : Float32[2.3787166f17, 9.364843f16, 3.5936449f16, 1.3888621f16, 5.4925405f15, 2.2399335f15]
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

```jldoctest
julia> T = collect(1073.0f0:40:1273.0f0);

julia> P = 2 .+ zero(T);

julia> dg = collect(3.0f0:8.0f-1:7.0f0);

julia> σ = collect(7.5f0:1.0f0:12.5f0) .* 1.0f-3;

julia> ϕ = collect(1.0f-2:2.0f-3:2.0f-2);

julia> T_solidus = 1200 + 273 .+ zero(T);

julia> model = xfit_premelt(T, P, dg, σ, ϕ, T_solidus)
Model : xfit_premelt
Temperature (K) : Float32[1073.0, 1113.0, 1153.0, 1193.0, 1233.0, 1273.0]
Pressure (GPa) : Float32[2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
grain size(μm) : Float32[3.0, 3.8, 4.6, 5.4, 6.2, 7.0]
Shear stress (GPa) : Float32[0.0075000003, 0.0085, 0.009500001, 0.010500001, 0.011500001, 0.0125]
Porosity : Float32[0.01, 0.012, 0.014, 0.016, 0.018, 0.02]
Solidus Temperature (K) : Float32[1473.0, 1473.0, 1473.0, 1473.0, 1473.0, 1473.0]

julia> forward(model, [])
Rock physics Response : RockphyViscous
Strain rate : Float32[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
Viscosity (Pa s) : Float32[7.6341115f18, 2.2587335f18, 6.6676584f17, 2.024383f17, 6.40973f16, 2.1291828f16]
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
