struct RockphyCondDistribution{T <: Union{Function, Nothing}} <:
       AbstractRockphyResponseDistribution
    σ::T
end

abstract type AbstractCondModelDistribution <: AbstractRockphyModelDistribution end

"""
    SEO3Distribution(T)

Model Distribution for `SEO3`[@ref].

## Arguments

  - `T` : Temperature of olivine (in K)

## Usage

Refer to the documentation for usage examples.
"""
mutable struct SEO3Distribution{F <: Union{Distribution, AbstractArray}} <:
               AbstractCondModelDistribution
    T::F
end

"""
    UHO2014Distribution(T, Ch2o_ol)

Model Distribution for `UHO2014`[@ref].

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_ol` : Water concentration in olivine (in ppm)

## Usage

Refer to the documentation for usage examples.
"""
mutable struct UHO2014Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray}} <:
               AbstractCondModelDistribution
    T::F1
    Ch2o_ol::F2
end

"""
    Jones2012Distribution(T, Ch2o_ol)

Model Distribution for `Jones2012`[@ref].

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_ol` : Water concentration in olivine (in ppm)

## Usage

Refer to the documentation for usage examples.
"""
mutable struct Jones2012Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray}} <:
               AbstractCondModelDistribution
    T::F1
    Ch2o_ol::F2
end

"""
    Poe2010Distribution(T, Ch2o_ol)

Model Distribution for `Poe2010`[@ref].

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_ol` : Water concentration in olivine (in ppm)

## Usage

Refer to the documentation for usage examples.
"""
mutable struct Poe2010Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray}} <:
               AbstractCondModelDistribution
    T::F1
    Ch2o_ol::F2
end

"""
    Wang2006Distribution(T, Ch2o_ol)

Model Distribution for `Wang2006`[@ref].

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_ol` : Water concentration in olivine (in ppm)

## Usage

Refer to the documentation for usage examples.
"""
mutable struct Wang2006Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray}} <:
               AbstractCondModelDistribution
    T::F1
    Ch2o_ol::F2
end

"""
    Yoshino2009Distribution(T, Ch2o_ol)

Model Distribution for `Yoshino2009`[@ref].

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_ol` : Water concentration in olivine (in ppm)

## Usage

Refer to the documentation for usage examples.
"""
mutable struct Yoshino2009Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray}} <:
               AbstractCondModelDistribution
    T::F1
    Ch2o_ol::F2
end

"""
    const_matrixDistribution(σ)

Model Distribution for `const_matrix`[@ref].

## Arguments

  - `σ` : Conductivity of the phase

## Usage

Refer to the documentation for usage examples.
"""
mutable struct const_matrixDistribution{F <: Union{Distribution, AbstractArray}} <:
               AbstractCondModelDistribution
    σ::F
end

"""
    Ni2011Distribution(T, Ch2o_m)

Model Distribution for `Ni2011`[@ref].

## Arguments

  - `T` : Temperature of melt (in K)
  - `Ch2o_m` : Water concentration in melt (in ppm)

## Usage

Refer to the documentation for usage examples.
"""
mutable struct Ni2011Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray}} <:
               AbstractCondModelDistribution
    T::F1
    Ch2o_m::F2
end

"""
    Sifre2014Distribution(T, Ch2o_m, Cco2_m)

Model Distribution for `Sifre2014`[@ref].

## Arguments

  - `T` : Temperature of melt (in K)
  - `Ch2o_m` : Water concentration in melt (in ppm)
  - `Cco2_m` : CO₂ concentration in melt (in ppm)

## Usage

Refer to the documentation for usage examples.
"""
mutable struct Sifre2014Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray},
    F3 <: Union{Distribution, AbstractArray}} <: AbstractCondModelDistribution
    T::F1
    Ch2o_m::F2
    Cco2_m::F3
end

"""
    Gaillard2008Distribution(T)

Model Distribution for `Gaillard2008`[@ref].

## Arguments

  - `T` : Temperature of melt (in K)

## Usage

Refer to the documentation for usage examples.
"""
mutable struct Gaillard2008Distribution{F <: Union{Distribution, AbstractArray}} <:
               AbstractCondModelDistribution
    T::F
end

"""
    Dai_Karato2009Distribution(T, Ch2o_opx)

Electrical conductivity model for olivine dependent on temperature and water concentration.

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_opx` : water concentration in orthopyroxene (in ppm)

## Usage

Refer to the documentation for usage examples.
"""
mutable struct Dai_Karato2009Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray}} <:
               AbstractCondModelDistribution
    T::F1
    Ch2o_opx::F2
end

"""
    Zhang2012Distribution(T, Ch2o_opx)

Electrical conductivity model for olivine dependent on temperature and water concentration.

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_opx` : water concentration in olivine (in ppm)

## Usage

Refer to the documentation for usage examples.
"""
mutable struct Zhang2012Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray}} <:
               AbstractCondModelDistribution
    T::F1
    Ch2o_opx::F2
end

# ========================================================================================================== 
# clinopyroxene

"""
    Yang2011Yang2011Distribution(T, Ch2o_cpx)

Electrical conductivity model for olivine dependent on temperature and water concentration.

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_cpx` : water concentration in olivine (in ppm)

## Usage

Refer to the documentation for usage examples.
"""
mutable struct Yang2011Distribution{
    F1 <: Union{Distribution, AbstractArray}, F2 <: Union{Distribution, AbstractArray}} <:
               AbstractCondModelDistribution
    T::F1
    Ch2o_cpx::F2
end
