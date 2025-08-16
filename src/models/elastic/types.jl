abstract type AbstractElasticModel <: AbstractRockphyModel end

## response

mutable struct RockphyElastic{T1, T2, T3, T4} <: AbstractRockphyResponse
    G::T1
    K::T2
    Vp::T3
    Vs::T4
end

"""
    anharmonic(T, P, ρ)

Calculate unrelaxed elastic bulk moduli, shear moduli, Vp and Vs
using anharmonic scaling

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `ρ` : Density of the rock (kg/m³)

## Keyword Arguments

    - `params` : Various coefficients required for calculation. 
    Available options are 
        - `params_anharmonic.Isaak1992`
        - `params_anharmonic.Cammarono2003`
    Defaults to the ones provided by Isaak (1992), check references. 
    To investigate coefficients, call `default_params(Val{anharmonic}())`. 
    To modify coefficients, check the relevant documentation page.

## Usage

```julia
T = collect(1273.0f0:30:1573.0f0)
P = 2 .+ zero(T)
ρ = collect(3300.0f0:100.0f0:4300.0f0)

model = anharmonic(T, P, ρ)

forward(model, [])
```

## References

  - Cammarano et al. (2003), "Inferring upper-mantle temperatures from seismic velocities",
    Physics of the Earth and Planetary Interiors, Volume 138, Issues 3–4,
    https://doi.org/10.1016/S0031-9201(03)00156-0

  - Isaak, D. G. (1992), "High‐temperature elasticity of iron‐bearing olivines",
    J. Geophys. Res., 97( B2), 1871– 1885,
    https://doi.org/10.1029/91JB02675
"""
mutable struct anharmonic{T1, T2, T3} <: AbstractElasticModel
    T::T1
    P::T2
    ρ::T3
end

"""
    anharmonic_poro(T, P, ρ, ϕ)

Calculate unrelaxed elastic bulk moduli, shear moduli, Vp and Vs
after applying correction for poro-elasticity using Takei (2002)
on anharmonic scaling

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `ρ` : Density of the rock (kg/m³)
    - `ϕ` : Porosity of the rock

## Keyword Arguments

    - `params` : Various coefficients required for calculation. 
    Contains a `p_anharmonic` which have the coefficients for [`anharmonic`](@ref)
    calculations. 
    To investigate coefficients, call `default_params(Val{anharmonic_poro}())`. 
    To modify coefficients, check the relevant documentation page.

## Usage

```julia
T = collect(1273.0f0:30:1573.0f0)
P = 2 .+ zero(T)
ρ = collect(3300.0f0:100.0f0:4300.0f0)
ϕ = collect(1.0f-2:1.0f-3:2.0f-2)

model = anharmonic_poro(T, P, ρ, ϕ)

forward(model, [])
```

## References

  - Takei, 2002, "Effect of pore geometry on VP/VS: From equilibrium geometry to crack",
    JGR Solid Earth,
    https://doi.org/10.1029/2001JB000522
"""
mutable struct anharmonic_poro{T1, T2, T3, T4} <: AbstractElasticModel
    T::T1
    P::T2
    ρ::T3
    ϕ::T4
end

"""
    SLB2005(T, P)

Calculate upper mantle shear velocity using
Stixrude and Lithgow‐Bertelloni (2005) fit of upper mantle Vs.

!!! warning


**Note that the other parameters (elastic moudli and Vp) are populated with zeros.**

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)

## Keyword Arguments

    - `params` : Various coefficients required for calculation. 
    Does not have any coefficients that can be changed and hence is an empty `NamedTuple`

## Usage

```julia
T = collect(1273.0f0:30:1573.0f0)
P = 2 .+ zero(T)

model = SLB2005(T, P)

forward(model, [])
```

## References

  - Stixrude and Lithgow‐Bertelloni (2005), "Mineralogy and elasticity of the oceanic upper mantle: Origin of the low‐velocity zone."
    JGR 110.B3,
    https://doi.org/10.1029/2004JB002965
"""
mutable struct SLB2005{T1, T2} <: AbstractElasticModel
    T::T1
    P::T2
end

default_params(::Type{T}) where {T <: anharmonic} = default_params_anharmonic
default_params(::Type{T}) where {T <: anharmonic_poro} = default_params_anharmonic_poro
default_params(::Type{T}) where {T <: SLB2005} = (;)
