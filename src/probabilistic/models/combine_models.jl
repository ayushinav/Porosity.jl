mutable struct multi_rp_responseDistribution{T1, T2, T3, T4} <:
               AbstractRockphyResponseDistribution
    cond::T1
    elastic::T2
    visc::T3
    anelastic::T4
end

"""
    multi_rp_modelDistributionType(con, elastic, visc, anelastic)

Rock physics model distribution type to capture multiple rock physics model distribution

## Arguments

  - `cond` : conductivity rock physics modeldistribution  type, subtype of `AbstractCondModelDistributionType`
  - `elastic` : elastic rock physics distribution model type, subtype of `AbstractElasticModelDistributionType`
  - `visc` : viscous rock physics distribution model type, subtype of `AbstractViscousModelDistributionType`
  - `anelastic` : anelastic rock physics distribution model type, subtype of `AbstractAnelasticModelDistributionType`

## Usage

Similar to `multi_rp_modelDistributionType` but instead accepts `Distribution`
or `Nothing`

```julia
multi_rp_modelDistributionType()(SEO3Distribution, anharmonicDistribution, HK2003Distribution, Nothing)
```

Pass `Nothing` for the types you do not want responses of, eg. above
does not compute for the `anelastic` type
"""
mutable struct multi_rp_modelDistributionType{T1, T2, T3, T4}
    cond::Type{T1}
    elastic::Type{T2}
    visc::Type{T3}
    anelastic::Type{T4}
end

function multi_rp_modelDistributionType(cond, elastic, visc, anelastic)
    cond_ = isa(cond, Type) ? cond : typeof(cond)
    elastic_ = isa(elastic, Type) ? elastic : typeof(elastic)
    visc_ = isa(visc, Type) ? visc : typeof(visc)
    anelastic_ = isa(anelastic, Type) ? anelastic : typeof(anelastic)
    return multi_rp_modelDistributionType(cond_, elastic_, visc_, anelastic_)
end

"""
    multi_rp_model(con, elastic, visc, anelastic)

Rock physics model to capture multiple rock physics model, susually constructed through `multi_rp_modelDistributionType`[@ref]

## Arguments

  - `cond` : conductivity rock physics modelDistribution
  - `elastic` : elastic rock physics modelDistribution
  - `visc` : viscous rock physics modelDistribution
  - `anelastic` : anelastic rock physics modelDistribution

## Usage

```julia
m = multi_rp_modelDistributionType()(SEO3Distribution, anharmonicDistribution, Nothing, Nothing)
ps_nt_dist = (; T=product_distribution(Uniform(1200.0f0, 1400.0f0)), P=[3.0f0],
    ρ=product_distribution(Uniform(80.0f0, 120.0f0)), ϕ=[0.1f0])
model = m(ps_nt_dist)
```
"""
mutable struct multi_rp_modelDistribution{T1, T2, T3, T4} <:
               AbstractRockphyModelDistribution
    cond::T1
    elastic::T2
    visc::T3
    anelastic::T4
end

function (model::multi_rp_modelDistributionType)(ps::NamedTuple)
    v = map((x) -> from_nt(getproperty(model, x), ps), propertynames(model))
    return multi_rp_modelDistribution(v...)
end
