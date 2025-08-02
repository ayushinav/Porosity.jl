"""
    two_phase_modelDistributionType(m1, m2, mix)

Rock physics model distribution type to combine two phases.

## Arguments

  - `m1` : model distribution type corresponding to phase 1
  - `m2` : model distribution type corresponding to phase 2
  - `mix` : mixing type, available options are `HS_1962_plus()`, `HS1962_minus`, `MAL(m)`

## Usage

```julia
two_phase_modelType(SEO3Distribution, Ni2011Distribution, HS1962_plus())
```
"""
mutable struct two_phase_modelDistributionType{T1, T2, M}
    m1::Type{T1}
    m2::Type{T2}
    mix::M
end

"""
    two_phase_modelDistribution(ϕ, m1, m2, mix)

Rock physics model distribution to combine two phases, usually constructed through `two_phase_modelDistributionType`[@ref]

## Arguments

  - `ϕ` : distribution (or value) of vol. fraction of the **second** phase
  - `m1` : model distribution corresponding to phase 1
  - `m2` : model distribution corresponding to phase 2
  - `mix` : mixing type, available options are `HS_1962_plus()`, `HS1962_minus`, `MAL(m)`

## Usage

```julia
m = two_phase_modelType(SEO3Distribution, Ni2011Distribution, HS1962_plus())
ps_nt_dist = (; T=product_distribution(Uniform(1200.0f0, 1400.0f0)),
    Ch2o_m=MvNormal([100.0f0], diagm([20.0f0])), ϕ=[0.1f0])
model = m(ps_nt_dist)

resp = forward(model)
```
"""
mutable struct two_phase_modelDistribution{V, T1, T2, M} <: AbstractRockphyModelDistribution
    ϕ::V
    m1::T1
    m2::T2
    mix::M
end

two_phase_modelDistributionType(m1) = m1
two_phase_modelDistributionType(m1, m::phase_mixing) = m1

function (model::two_phase_modelDistributionType)(ps::NamedTuple)
    mix = model.mix
    ϕ = ps.ϕ

    v1 = from_nt(model.m1, ps)
    v2 = from_nt(model.m2, ps)

    return two_phase_modelDistribution(ϕ, v1, v2, mix)
end

function from_nt(m::Type{T}, nt::NamedTuple) where {T <: two_phase_modelDistributionType}
    ϕ = nt.ϕ

    m1 = m.types[1].parameters[1]
    m2 = m.types[2].parameters[1]
    # mix = m.types[3] #.parameters[1]

    model1 = MT.from_nt(m1, nt)
    model2 = MT.from_nt(m2, nt)

    mix = from_nt(m.types[3], nt)

    return two_phase_modelDistribution(ϕ, model1, model2, mix)
end

# =========

mutable struct multi_phase_modelDistributionType{T1, T2, T3, T4, T5, T6, T7, T8, M}
    m1::Type{T1}
    m2::Type{T2}
    m3::Type{T3}
    m4::Type{T4}
    m5::Type{T5}
    m6::Type{T6}
    m7::Type{T7}
    m8::Type{T8}
    mix::M
end

multi_phase_modelDistributionType(m1) = m1
multi_phase_modelDistributionType(m1, m::phase_mixing) = m1

for i in 2:7
    args = [Symbol("m$k") for k in 1:i]
    last_args = :(m::phase_mixing)
    expr_lhs = Expr(:call, :multi_phase_modelDistributionType, args..., last_args)

    args2 = [Nothing for k in (i + 1):8]
    expr_rhs = Expr(:call, :multi_phase_modelDistributionType, args..., args2..., last_args)

    expr = Expr(:function, expr_lhs, expr_rhs)
    eval(expr)
end

mutable struct multi_phase_modelDistribution{V, T1, T2, T3, T4, T5, T6, T7, T8, M} <:
               AbstractRockphyModelDistribution
    ϕ::V
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

function (model::multi_phase_modelDistributionType)(ps::NamedTuple)
    pnames = propertynames(model)
    mix = model.mix
    ϕ = ps.ϕ

    # ϕ_vec = rearrange_ϕ(ps.ϕ, model)

    v1 = from_nt(getproperty(model, pnames[1]), ps)
    v2 = from_nt(getproperty(model, pnames[2]), ps)
    v3 = from_nt(getproperty(model, pnames[3]), ps)
    v4 = from_nt(getproperty(model, pnames[4]), ps)
    v5 = from_nt(getproperty(model, pnames[5]), ps)
    v6 = from_nt(getproperty(model, pnames[6]), ps)
    v7 = from_nt(getproperty(model, pnames[7]), ps)
    v8 = from_nt(getproperty(model, pnames[8]), ps)
    return multi_phase_modelDistribution(ϕ, v1, v2, v3, v4, v5, v6, v7, v8, mix)
end

function from_nt(m::Type{T}, nt::NamedTuple) where {T <: multi_phase_modelDistributionType}
    ϕ = nt.ϕ
    m1 = m.types[1].parameters[1]
    m2 = m.types[2].parameters[1]
    m3 = m.types[3].parameters[1]
    m4 = m.types[4].parameters[1]
    m5 = m.types[5].parameters[1]
    m6 = m.types[6].parameters[1]
    m7 = m.types[7].parameters[1]
    m8 = m.types[8].parameters[1]

    model1 = MT.from_nt(m1, nt)
    model2 = MT.from_nt(m2, nt)
    model3 = MT.from_nt(m3, nt)
    model4 = MT.from_nt(m4, nt)
    model5 = MT.from_nt(m5, nt)
    model6 = MT.from_nt(m6, nt)
    model7 = MT.from_nt(m7, nt)
    model8 = MT.from_nt(m8, nt)

    mix = from_nt(m.types[9], nt)

    return multi_phase_modelDistribution(
        ϕ, model1, model2, model3, model4, model5, model6, model7, model8, mix)
end
