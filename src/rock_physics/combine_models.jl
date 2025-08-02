mutable struct multi_rp_response{T1, T2, T3, T4} <: AbstractRockphyResponse
    cond::T1
    elastic::T2
    visc::T3
    anelastic::T4
end

"""
    multi_rp_modelType(con, elastic, visc, anelastic)

Rock physics model type to capture multiple rock physics models

## Arguments

  - `cond` : conductivity rock physics model type, subtype of `AbstractCondModel`
  - `elastic` : elastic rock physics model type, subtype of `AbstractElasticModel`
  - `visc` : viscous rock physics model type, subtype of `AbstractViscousModel`
  - `anelastic` : anelastic rock physics model type, subtype of `AbstractAnelasticModel`

## Usage

```julia
multi_rp_modelType()(SEO3, anharmonic, HK2003, Nothing)
```

Pass `Nothing` for the types you do not want responses of, eg. above
does not compute for the `anelastic` type
"""
mutable struct multi_rp_modelType{T1, T2, T3, T4}
    cond::Type{T1}
    elastic::Type{T2}
    visc::Type{T3}
    anelastic::Type{T4}
end

function multi_rp_modelType(cond, elastic, visc, anelastic)
    cond_ = isa(cond, Type) ? cond : typeof(cond)
    elastic_ = isa(elastic, Type) ? elastic : typeof(elastic)
    visc_ = isa(visc, Type) ? visc : typeof(visc)
    anelastic_ = isa(anelastic, Type) ? anelastic : typeof(anelastic)
    return multi_rp_modelType(cond_, elastic_, visc_, anelastic_)
end

"""
    multi_rp_model(cond, elastic, visc, anelastic)

Rock physics model to capture multiple rock physics models,
usually constructed through `multi_rp_modelType`[@ref]

## Arguments

  - `cond` : conductivity rock physics model
  - `elastic` : elastic rock physics model
  - `visc` : viscous rock physics model
  - `anelastic` : anelastic rock physics model

## Usage

```julia
m = multi_rp_modelType()(SEO3, anharmonic, Nothing, Nothing)
ps_nt = ps_nt = (;
    T=[800.0f0, 1000.0f0] .+ 273, P=3.0f0, ρ=3300.0f0, Ch2o_m=1000.0f0, ϕ=0.1f0)
model = m(ps_nt)

resp = forward(model)
```
"""
mutable struct multi_rp_model{T1, T2, T3, T4} <: AbstractRockphyModel
    cond::T1
    elastic::T2
    visc::T3
    anelastic::T4
end

function (model::multi_rp_modelType)(ps::NamedTuple)
    v = map((x) -> from_nt(getproperty(model, x), ps), propertynames(model))
    return multi_rp_model(v...)
end

function forward(model::multi_rp_model{T1, T2, T3, T4}, p, params) where {T1, T2, T3, T4}
    resp_cond = forward(model.cond, [], params.cond)
    resp_elastic = forward(model.elastic, [], params.elastic)
    resp_visc = forward(model.visc, [], params.visc)
    resp_anelastic = forward(model.anelastic, [], params.anelastic)

    multi_rp_response(resp_cond, resp_elastic, resp_visc, resp_anelastic)
end

function forward(model::multi_rp_model{T1, T2, T3, T4}, p) where {T1, T2, T3, T4}
    resp_cond = forward(model.cond, [])
    resp_elastic = forward(model.elastic, [])
    resp_visc = forward(model.visc, [])
    resp_anelastic = forward(model.anelastic, [])

    multi_rp_response(resp_cond, resp_elastic, resp_visc, resp_anelastic)
end

function default_params(::Type{multi_rp_modelType{T1, T2, T3, T4}}) where {T1, T2, T3, T4}
    (;
        zip([:cond, :elastic, :visc, :anelastic],
            [default_params(T1), default_params(T2),
                default_params(T3), default_params(T4)])...)
end

function from_nt(m::Type{T}, nt::NamedTuple) where {T <: multi_rp_modelType}
    m1 = T.parameters[1]
    m2 = T.parameters[2]
    m3 = T.parameters[3]
    m4 = T.parameters[4]

    model1 = MT.from_nt(m1, nt)
    model2 = MT.from_nt(m2, nt)
    model3 = MT.from_nt(m3, nt)
    model4 = MT.from_nt(m4, nt)

    return multi_rp_model(model1, model2, model3, model4)
end
