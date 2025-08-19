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

```jldoctest
julia> T = collect(1273.0f0:40:1473.0f0);

julia> P = 2 .+ zero(T);

julia> ρ = collect(3300.0f0:200.0f0:4300.0f0);

julia> model = anharmonic(T, P, ρ)
Model : anharmonic
Temperature (K) : Float32[1273.0, 1313.0, 1353.0, 1393.0, 1433.0, 1473.0]
Pressure (GPa) : Float32[2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
Density (kg/m³) : Float32[3300.0, 3500.0, 3700.0, 3900.0, 4100.0, 4300.0]

julia> forward(model, [])
Rock physics Response : RockphyElastic
Elastic shear modulus (Pa) : Float32[7.136702f10, 7.082302f10, 7.027902f10, 6.9735014f10, 6.919102f10, 6.864702f10]
Elastic bulk modulus (Pa) : Float32[1.1894502f11, 1.18038364f11, 1.1713171f11, 1.1622502f11, 1.1531837f11, 1.1441169f11]
Elastic P-wave velocity (m/s) : Float32[8054.7563, 7791.3696, 7548.708, 7324.0913, 7115.3057, 6920.496]
Elastic S-wave velocity (m/s) : Float32[4650.416, 4498.3496, 4358.2485, 4228.5664, 4108.0234, 3995.5503]
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

```jldoctest
julia> T = collect(1273.0f0:40:1473.0f0);

julia> P = 2 .+ zero(T);

julia> ρ = collect(3300.0f0:200.0f0:4300.0f0);

julia> ϕ = collect(1.0f-2:2.0f-3:2.0f-2);

julia> model = anharmonic_poro(T, P, ρ, ϕ)
Model : anharmonic_poro
Temperature (K) : Float32[1273.0, 1313.0, 1353.0, 1393.0, 1433.0, 1473.0]
Pressure (GPa) : Float32[2.0, 2.0, 2.0, 2.0, 2.0, 2.0]
Density (kg/m³) : Float32[3300.0, 3500.0, 3700.0, 3900.0, 4100.0, 4300.0]
Porosity : Float32[0.01, 0.012, 0.014, 0.016, 0.018, 0.02]

julia> forward(model, [])
Rock physics Response : RockphyElastic
Elastic shear modulus (Pa) : Float32[6.90533f10, 6.8138574f10, 6.723185f10, 6.633203f10, 6.5438487f10, 6.4550883f10]
Elastic bulk modulus (Pa) : Float32[1.1665774f11, 1.1534928f11, 1.1405562f11, 1.127763f11, 1.11511036f11, 1.102595f11]
Elastic P-wave velocity (m/s) : Float32[7953.0596, 7675.5776, 7419.807, 7182.9395, 6962.659, 6757.0347]
Elastic S-wave velocity (m/s) : Float32[4574.4116, 4412.2744, 4262.7188, 4124.1016, 3995.0728, 3874.5107]
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

```jldoctest
julia> T = collect(1273.0f0:40:1473.0f0);

julia> P = 2 .+ zero(T);

julia> model = SLB2005(T, P)
Model : SLB2005
Temperature (K) : Float32[1273.0, 1313.0, 1353.0, 1393.0, 1433.0, 1473.0]
Pressure (GPa) : Float32[2.0, 2.0, 2.0, 2.0, 2.0, 2.0]

julia> forward(model, [])
Rock physics Response : RockphyElastic
Elastic shear modulus (Pa) : Float32[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
Elastic bulk modulus (Pa) : Float32[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
Elastic P-wave velocity (m/s) : Float32[0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
Elastic S-wave velocity (m/s) : [4478.205990493298, 4463.085990846157, 4447.965991199017, 4432.845991551876, 4417.725991904736, 4402.605992257595]
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
