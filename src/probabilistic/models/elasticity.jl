struct RockphyElasticDistribution{
    T1 <: Union{Function, Nothing}, T2 <: Union{Function, Nothing},
    T3 <: Union{Function, Nothing}, T4 <: Union{Function, Nothing}} <:
       AbstractRockphyResponseDistribution
    G::T1
    K::T2
    Vp::T3
    Vs::T4
end

abstract type AbstractElasticModelDistribution <: AbstractRockphyModelDistribution end

"""
    anharmonicDistribution(T, P, ρ)

Model Distribution for `anharmonic_poro`[@ref].

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `ρ` : Density of the rock (kg/m³)

## Usage

Refer to the documentation for usage examples.
"""
mutable struct anharmonicDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}} <: AbstractElasticModelDistribution
    T::T1
    P::T2
    ρ::T3
end

"""
    anharmonic_poroDistribution(T, P, ρ, ϕ)

Model Distribution for `anharmonic`[@ref].

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `ρ` : Density of the rock (kg/m³)
    - `ϕ` : Porosity of the rock

## Usage

Refer to the documentation for usage examples.
"""
mutable struct anharmonic_poroDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray}} <:
               AbstractElasticModelDistribution
    T::T1
    P::T2
    ρ::T3
    ϕ::T4
end

"""
    SLB2005Distribution(T, P)

Model Distribution for `SLB2005`[@ref].

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)

## Usage

Refer to the documentation for usage examples.
"""
mutable struct SLB2005Distribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray}} <:
               AbstractElasticModelDistribution
    T::T1
    P::T2
end
