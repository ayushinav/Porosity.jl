"""
    two_phase_modelType(m1, m2, mix)

Rock physics model Type to combine two phases.

## Arguments

  - `m1` : model type corresponding to phase 1
  - `m2` : model type corresponding to phase 2
  - `mix` : mixing type, available options are `HS_1962_plus`, `HS1962_minus`, `MAL`

## Usage

```jldoctest; output = false
two_phase_modelType(SEO3, Ni2011, HS1962_plus)

# output

two_phase_modelType{SEO3, Ni2011, HS1962_plus}(SEO3, Ni2011, HS1962_plus)
```
"""
mutable struct two_phase_modelType{T1, T2, M}
    m1::Type{T1}
    m2::Type{T2}
    mix::Type{M}
end

"""
    two_phase_model(ϕ, m1, m2, mix)

Rock physics model to combine two phases, usually constructed through `two_phase_modelType`[@ref]

## Arguments

  - `ϕ` : Vol. fraction of the **second** phase
  - `m1` : model corresponding to phase 1
  - `m2` : model corresponding to phase 2
  - `mix` : mixing type, available options are `HS_1962_plus()`, `HS1962_minus`, `MAL(m)`

## Usage

```jldoctest
julia> m = two_phase_modelType(SEO3, Ni2011, HS1962_plus);

julia> T = [900.0f0, 1000.0f0] .+ 273;

julia> P = 3.0f0;

julia> ρ = 3300.0f0;

julia> Ch2o_m = 1000.0f0;

julia> ϕ = 0.1f0;

julia> ps_nt = ps_nt = (; T, P, ρ, Ch2o_m, ϕ);

julia> model = m(ps_nt)
Two phase composition using HS1962_plus()

*    m₁ (solid phase) : 
Model : SEO3
Temperature (K) : Float32[1173.0, 1273.0]
ϕ : 0.9

*    m₂ (liquid phase) : 
Model : SEO3
Temperature (K) : Float32[1173.0, 1273.0]
ϕ : 0.1

julia> resp = forward(model, [])
Rock physics Response : RockphyCond
log₁₀ conductivity (S/m) : Float32[-5.800449, -4.128538]
```
"""
mutable struct two_phase_model{V, T1, T2, M} <: AbstractRockphyModel
    ϕ::V
    m1::T1
    m2::T2
    mix::M
end

two_phase_modelType(m1) = m1
two_phase_modelType(m1, mix::Type{T}) where {T <: phase_mixing} = m1

function (model::two_phase_modelType)(ps::NamedTuple)
    ϕ = ps.ϕ
    v1 = from_nt(getproperty(model, :m1), ps)
    v2 = from_nt(getproperty(model, :m2), ps)
    mix = from_nt(model.mix, ps)
    return two_phase_model(ϕ, v1, v2, mix)
end

function SubsurfaceCore.forward(model::two_phase_model{V, T1, T2, M},
        p) where {
        V, M <: two_phase_mix_types, T1 <: AbstractCondModel, T2 <: AbstractCondModel}
    σ1 = forward(model.m1, []).σ
    σ2 = forward(model.m2, []).σ

    @. σ1 = exp10(σ1)
    @. σ2 = exp10(σ2)

    σ = broadcast(
        (sig1, sig2, phi) -> mix_models([sig1, sig2], phi, model.mix), σ1, σ2, model.ϕ)

    return RockphyCond(log10.(σ))
end

function SubsurfaceCore.forward(model::two_phase_model{V, T1, T2, M}, p,
        params) where {V, M, T1 <: AbstractCondModel, T2 <: AbstractCondModel}
    σ1 = forward(model.m1, [], params.m1).σ
    σ2 = forward(model.m2, [], params.m2).σ

    @. σ1 = exp10(σ1)
    @. σ2 = exp10(σ2)

    σ = broadcast(
        (sig1, sig2, phi) -> mix_models([sig1, sig2], phi, model.mix), σ1, σ2, model.ϕ)

    return RockphyCond(log10.(σ))
end

function SubsurfaceCore.default_params(::Type{two_phase_modelType{
        T1, T2, M}}) where {T1, T2, M}
    (; zip([:m1, :m2], [default_params(T1), default_params(T2)])...)
end

function SubsurfaceCore.default_params(::Type{two_phase_model{
        V, T1, T2, M}}) where {V, T1, T2, M}
    (; zip([:m1, :m2], [default_params(T1), default_params(T2)])...)
end

# following is needed for combine_models

# for constructing mdist
function SubsurfaceCore.from_nt(m::Type{T}, nt::NamedTuple) where {T <: two_phase_modelType}
    ϕ = nt.ϕ
    m1 = m.types[1].parameters[1]
    m2 = m.types[2].parameters[1]

    model1 = from_nt(m1, nt)
    model2 = from_nt(m2, nt)
    # @show m.types[3]
    mix = from_nt(m.types[3].parameters[1], nt)

    return two_phase_model(ϕ, model1, model2, mix)
end

# ==============================================================================
# multi-phase 

"""
    multi_phase_modelType(m1, m2, m3, m4, m5, m6, m7, m8, mix)
    multi_phase_modelType(m1, ..., mix)

Rock physics model Type to combine multiple (upto 8) phases.

## Arguments

  - `m1` : model type corresponding to phase 1
  - `m2` : model type corresponding to phase 2 (optional)
  - `m3` : model type corresponding to phase 3 (optional)
  - `m4` : model type corresponding to phase 4 (optional)
  - `m5` : model type corresponding to phase 5 (optional)
  - `m6` : model type corresponding to phase 6 (optional)
  - `m7` : model type corresponding to phase 7 (optional)
  - `m8` : model type corresponding to phase 8 (optional)
  - `mix` : mixing type, available options are `HS_plus_multi_phase`, `HS_minus_multi_phase`, `GAL`

## Usage

```julia
multi_phase_modelType(SEO3, Ni2011, Yang2011, HS_plus_multi_phase)
```
"""
mutable struct multi_phase_modelType{T1, T2, T3, T4, T5, T6, T7, T8, M}
    m1::Type{T1}
    m2::Type{T2}
    m3::Type{T3}
    m4::Type{T4}
    m5::Type{T5}
    m6::Type{T6}
    m7::Type{T7}
    m8::Type{T8}
    mix::Type{M}
end

multi_phase_modelType(m1) = m1
multi_phase_modelType(m1, m::Type{phase_mixing}) = m1

for i in 2:7
    args = [Symbol("m$k") for k in 1:i]
    T_ = :T
    last_args = :(m::Type{$T_})
    expr_call = Expr(:call, :multi_phase_modelType, args..., last_args)
    expr_where = :($T_ <: phase_mixing)
    expr_lhs = Expr(:where, expr_call, expr_where)

    args2 = [Nothing for k in (i + 1):8]
    expr_rhs = Expr(:call, :multi_phase_modelType, args..., args2..., :m)

    expr = Expr(:function, expr_lhs, expr_rhs)
    eval(expr)
end

"""
    multi_phase_model(m1, m2, m3, m4, m5, m6, m7, m8, mix)
    multi_phase_model(m1, ..., mix)

Rock physics model Type to combine multiple (upto 8) phases.

## Arguments

  - `m1` : model corresponding to phase 1
  - `m2` : model corresponding to phase 2 (optional)
  - `m3` : model corresponding to phase 3 (optional)
  - `m4` : model corresponding to phase 4 (optional)
  - `m5` : model corresponding to phase 5 (optional)
  - `m6` : model corresponding to phase 6 (optional)
  - `m7` : model corresponding to phase 7 (optional)
  - `m8` : model corresponding to phase 8 (optional)
  - `mix` : mixing type, available options are `HS_plus_multi_phase()`, `HS_minus_multi_phase`, `GAL(m)`

## Usage

```jldoctest
julia> m = multi_phase_modelType(SEO3, Ni2011, Yang2011, HS_plus_multi_phase);

julia> T = [900.0f0, 1000.0f0] .+ 273;

julia> P = 3.0f0;

julia> ρ = 3300.0f0;

julia> Ch2o_m = 1000.0f0;

julia> Ch2o_cpx = 100.0f0;

julia> ϕ = [0.8f0, 0.1f0];

julia> ps_nt = ps_nt = (; T, P, ρ, Ch2o_m, Ch2o_cpx, ϕ);

julia> model = m(ps_nt)
multi-phase composition using HS_plus_multi_phase()

*    m₁ : 
Model : SEO3
Temperature (K) : Float32[1173.0, 1273.0]
ϕ : 0.8

*    m₂ : 
Model : Ni2011
Temperature (K) : Float32[1173.0, 1273.0]
Water concentration in melt (ppm) : 1000.0
ϕ : 0.1

*    m₃ : 
Model : Yang2011
Temperature (K) : Float32[1173.0, 1273.0]
Water concentration in clinopyroxene (ppm) : 100.0
ϕ : 0.099999964

julia> forward(model, [])
Rock physics Response : RockphyCond
log₁₀ conductivity (S/m) : Float32[-26.768837, -3.9526708]
```
"""
mutable struct multi_phase_model{T, T1, T2, T3, T4, T5, T6, T7, T8, M} <:
               AbstractRockphyModel
    ϕ::T
    m1::T1
    m2::T2
    m3::T3
    m4::T4
    m5::T5
    m6::T6
    m7::T7
    m8::T8
    mix::M
end

function rearrange_ϕ(ϕ, model::multi_phase_modelType)
    @assert sum(ϕ)≤1 "Σϕᵢ = $(sum(ϕ)) should be ≤ 1."

    fnames = propertynames(model)[1:(end - 1)]
    fnames = filter(f -> getfield(model, f) !== Nothing, fnames)
    msg = """
    Vol. frac of last component is defined automatically once the others are defined. 
    Make sure the length of porosity parameter `ϕ` is one less than the length of number of components.
    length of ϕ = $(length(ϕ))
    length of components = $(length(fnames))
    Check out the relevant documentation.
    """

    @assert length(ϕ)==(length(fnames) - 1) msg

    c = length(ϕ)
    ϕ_vec = zeros(eltype(ϕ), length(fnames))

    ϕ_vec[1:c] .= ϕ
    ϕ_vec[c + 1] = 1 - sum(ϕ)
    return ϕ_vec
end

function (model::multi_phase_modelType)(ps::NamedTuple)
    mix = from_nt(model.mix, ps)

    ϕ_vec = rearrange_ϕ(ps.ϕ, model)

    v1 = from_nt(model.m1, ps)
    v2 = from_nt(model.m2, ps)
    v3 = from_nt(model.m3, ps)
    v4 = from_nt(model.m4, ps)
    v5 = from_nt(model.m5, ps)
    v6 = from_nt(model.m6, ps)
    v7 = from_nt(model.m7, ps)
    v8 = from_nt(model.m8, ps)
    return multi_phase_model(ϕ_vec, v1, v2, v3, v4, v5, v6, v7, v8, mix)
end

function SubsurfaceCore.forward(
        model::multi_phase_model{V, T1, T2, T3, T4, T5, T6, T7, T8, M},
        p) where {V, M, T1 <: AbstractCondModel, T2, T3, T4, T5, T6, T7, T8}
    σ1 = (isnothing(model.m1)) ? nothing : forward(model.m1, []).σ .|> exp10
    σ2 = (isnothing(model.m2)) ? nothing : forward(model.m2, []).σ .|> exp10
    σ3 = (isnothing(model.m3)) ? nothing : forward(model.m3, []).σ .|> exp10
    σ4 = (isnothing(model.m4)) ? nothing : forward(model.m4, []).σ .|> exp10
    σ5 = (isnothing(model.m5)) ? nothing : forward(model.m5, []).σ .|> exp10
    σ6 = (isnothing(model.m6)) ? nothing : forward(model.m6, []).σ .|> exp10
    σ7 = (isnothing(model.m7)) ? nothing : forward(model.m7, []).σ .|> exp10
    σ8 = (isnothing(model.m8)) ? nothing : forward(model.m8, []).σ .|> exp10

    σ_tup_ = (σ1, σ2, σ3, σ4, σ5, σ6, σ7, σ8)
    σ_tup = filter(f -> !isnothing(f), σ_tup_) |> Tuple

    σ = broadcast_helper_(σ_tup, model.ϕ, model.mix, Val{length(σ_tup)}())
    return RockphyCond(log10.(σ))
end

function SubsurfaceCore.forward(
        model::multi_phase_model{V, T1, T2, T3, T4, T5, T6, T7, T8, M}, p,
        params) where {V, M, T1 <: AbstractCondModel, T2, T3, T4, T5, T6, T7, T8}
    σ1 = (isnothing(model.m1)) ? Nothing : forward(model.m1, params.m1).σ .|> exp10
    σ2 = (isnothing(model.m2)) ? Nothing : forward(model.m2, params.m2).σ .|> exp10
    σ3 = (isnothing(model.m3)) ? Nothing : forward(model.m3, params.m3).σ .|> exp10
    σ4 = (isnothing(model.m4)) ? Nothing : forward(model.m4, params.m4).σ .|> exp10
    σ5 = (isnothing(model.m5)) ? Nothing : forward(model.m5, params.m5).σ .|> exp10
    σ6 = (isnothing(model.m6)) ? Nothing : forward(model.m6, params.m6).σ .|> exp10
    σ7 = (isnothing(model.m7)) ? Nothing : forward(model.m7, params.m7).σ .|> exp10
    σ8 = (isnothing(model.m8)) ? Nothing : forward(model.m8, params.m8).σ .|> exp10

    σ_tup_ = (σ1, σ2, σ3, σ4, σ5, σ6, σ7, σ8)
    σ_tup = filter(f -> isa(f, AbstractArray), σ_tup_) |> Tuple

    σ = broadcast_helper_(σ_tup, model.ϕ, model.mix, Val{length(σ_tup)}())

    return RockphyCond(log10.(σ))
end

function broadcast_helper_(m_tup, phi, mix, ::Val{2})
    broadcast((m1, m2) -> mix_models((m1, m2), phi, mix), m_tup[1], m_tup[2])
end

function broadcast_helper_(m_tup, phi, mix, ::Val{3})
    broadcast(
        (m1, m2, m3) -> mix_models((m1, m2, m3), phi, mix), m_tup[1], m_tup[2], m_tup[3])
end

function broadcast_helper_(m_tup, phi, mix, ::Val{4})
    broadcast((m1, m2, m3, m4) -> mix_models((m1, m2, m3, m4), phi, mix),
        m_tup[1], m_tup[2], m_tup[3], m_tup[4])
end

function broadcast_helper_(m_tup, phi, mix, ::Val{5})
    broadcast((m1, m2, m3, m4, m5) -> mix_models((m1, m2, m3, m4, m5), phi, mix),
        m_tup[1], m_tup[2], m_tup[3], m_tup[4], m_tup[5])
end

function broadcast_helper_(m_tup, phi, mix, ::Val{6})
    broadcast((m1, m2, m3, m4, m5, m6) -> mix_models((m1, m2, m3, m4, m5, m6), phi, mix),
        m_tup[1], m_tup[2], m_tup[3], m_tup[4], m_tup[5], m_tup[6])
end

function broadcast_helper_(m_tup, phi, mix, ::Val{7})
    broadcast(
        (m1, m2, m3, m4, m5, m6, m7) -> mix_models((m1, m2, m3, m4, m5, m6, m7), phi, mix),
        m_tup[1], m_tup[2], m_tup[3], m_tup[4], m_tup[5], m_tup[6], m_tup[7])
end

function broadcast_helper_(m_tup, phi, mix, ::Val{8})
    broadcast(
        (m1, m2, m3, m4, m5, m6, m7, m8) -> mix_models(
            (m1, m2, m3, m4, m5, m6, m7, m8), phi, mix),
        m_tup[1],
        m_tup[2],
        m_tup[3],
        m_tup[4],
        m_tup[5],
        m_tup[6],
        m_tup[7],
        m_tup[8])
end

function default_params(::Type{multi_phase_modelType{
        T1, T2, T3, T4, T5, T6, T7, T8, M}}) where {T1, T2, T3, T4, T5, T6, T7, T8, M}
    (;
        zip([:m1, :m2, :m3, :m4, :m5, :m6, :m7, :m8],
            [default_params(T1), default_params(T2), default_params(T3),
                default_params(T4), default_params(T5), default_params(T6),
                default_params(T7), default_params(T8)])...)
end

function SubsurfaceCore.from_nt(
        m::Type{T}, nt::NamedTuple) where {T <: multi_phase_modelType}
    ϕ = nt.ϕ
    m1 = m.types[1].parameters[1]
    m2 = m.types[2].parameters[1]
    m3 = m.types[3].parameters[1]
    m4 = m.types[4].parameters[1]
    m5 = m.types[5].parameters[1]
    m6 = m.types[6].parameters[1]
    m7 = m.types[7].parameters[1]
    m8 = m.types[8].parameters[1]

    model1 = from_nt(m1, nt)
    model2 = from_nt(m2, nt)
    model3 = from_nt(m3, nt)
    model4 = from_nt(m4, nt)
    model5 = from_nt(m5, nt)
    model6 = from_nt(m6, nt)
    model7 = from_nt(m7, nt)
    model8 = from_nt(m8, nt)

    mix = from_nt(m.types[9].parameters[1], nt)

    ϕ_vec = rearrange_ϕ(
        ϕ, multi_phase_modelType(m1, m2, m3, m4, m5, m6, m7, m8, m.types[9].parameters[1]))

    return multi_phase_model(
        ϕ_vec, model1, model2, model3, model4, model5, model6, model7, model8, mix)
end
