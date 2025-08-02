struct RockphyViscousDistribution{
    T1 <: Union{Function, Nothing}, T2 <: Union{Function, Nothing}} <:
       AbstractRockphyResponseDistribution
    ϵ_rate::T1
    η::T2
end

"""
    HZK2011Distribution(T, P, dg, σ, ϕ)

Model Distribution for `HZK2011`[@ref].

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `dg`: Grain size (μm)
    - `σ` : Shear stress (GPa)
    - `ϕ` : Porosity

## Usage

Refer to the documentation for usage examples.
"""
mutable struct HZK2011Distribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}} <: AbstractRockphyModelDistribution
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
end

"""
    HK2003(T, P, dg, σ, ϕ, Ch2o_ol)

Model Distribution for `HK2003`[@ref].

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `dg`: Grain size (μm)
    - `σ` : Shear stress (GPa)
    - `ϕ` : Porosity
    - `Ch2o_ol` : Water concentration in olivine (ppm), defaults to 0 ppm.

## Usage

Refer to the documentation for usage examples.
"""
mutable struct HK2003Distribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}, T6 <: Union{Distribution, AbstractArray}} <:
               AbstractRockphyModelDistribution
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    Ch2o_ol::T6
end

"""
    xfit_premeltDistribution(T, P, dg, σ, ϕ, Ch2o_ol)

Model Distribution for `xfit_premelt[@ref].

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `dg`: Grain size (μm)
    - `σ` : Shear stress (GPa)
    - `ϕ` : Porosity
    - `T_solidus` : Solidus temperature

## Usage

Refer to the documentation for usage examples.
"""
mutable struct xfit_premeltDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray},
    T3 <: Union{Distribution, AbstractArray}, T4 <: Union{Distribution, AbstractArray},
    T5 <: Union{Distribution, AbstractArray}, T6 <: Union{Distribution, AbstractArray}} <:
               AbstractRockphyModelDistribution
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    T_solidus::T6
end
