# utils to help bump for Abstract types

import LinearAlgebra: zero

function zero(x::resp) where {resp <: AbstractResponse}
    typeof(x)([zero(getfield(x, k)) for k in fieldnames(resp)]...)
end
function zero(x::model) where {model <: AbstractModel}
    typeof(x)([zero(getfield(x, k)) for k in fieldnames(model)]...)
end

import Base: copy
function copy(x::resp) where {resp <: AbstractResponse}
    typeof(x)([copy(getfield(x, k)) for k in fieldnames(resp)]...)
end
function copy(x::model) where {model <: AbstractModel}
    typeof(x)([copy(getfield(x, k)) for k in fieldnames(model)]...)
end

# NamedTuple manipulation

@generated function from_nt(::Type{T}, nt::NamedTuple) where {T}
    # @show T
    fnames = fieldnames(T)
    args = [:(getproperty(nt, $(QuoteNode(f)))) for f in fnames]
    return :(T($(args...)))
end

function from_nt(::Type{Nothing}; nt::NamedTuple)
    (;)
end

function to_nt(s)
    T = typeof(s)
    names = fieldnames(T)
    vals = ntuple(i -> getfield(s, names[i]), length(names))
    NamedTuple{names}(vals)
end

to_resp_nt(d::T) where {T <: AbstractResponse} = to_nt(d)

function to_resp_nt(d::T) where {T <: multi_rp_response}
    return merge(to_nt(d.cond), to_nt(d.elastic), to_nt(d.visc), to_nt(d.anelastic))
end

# forward manipulation

forward(m::Nothing, p, params=(;)) = nothing
default_params(::Type{Nothing}) = (;)

function forward_helper(
        m::Type{T}, m0, vars, response_trans_utils, params) where {T <: AbstractGeophyModel}
    model = from_nt(m, m0)
    resp_nt = to_resp_nt(forward(model, vars, response_trans_utils))
    return resp_nt
end

function forward_helper(m::Type{T}, m0, vars, response_trans_utils,
        params) where {T <: AbstractRockphyModel}
    model = from_nt(m, m0)
    resp_nt = to_resp_nt(forward(model, vars, params))
    for k in propertynames(resp_nt)
        broadcast!(response_trans_utils[k].tf, resp_nt[k], resp_nt[k])
    end
    return resp_nt
end

function forward_helper(
        m::Type{T}, m0, vars, response_trans_utils, params) where {T <: two_phase_modelType}
    model = from_nt(m, m0)
    resp_nt = to_resp_nt(forward(model, vars, params))
    for k in propertynames(resp_nt)
        broadcast!(response_trans_utils[k].tf, resp_nt[k], resp_nt[k])
    end
    return resp_nt
end

function forward_helper(m::Type{T}, m0, vars, response_trans_utils,
        params) where {T <: multi_phase_modelType}
    model = from_nt(m, m0)
    resp_nt = to_resp_nt(forward(model, vars, params))
    for k in propertynames(resp_nt)
        broadcast!(response_trans_utils[k].tf, resp_nt[k], resp_nt[k])
    end
    return resp_nt
end

function forward_helper(
        m::Type{T}, m0, vars, response_trans_utils, params) where {T <: multi_rp_modelType}
    model = from_nt(m, m0)
    resp_nt = to_resp_nt(forward(model, vars, params))
    for k in propertynames(resp_nt)
        broadcast!(response_trans_utils[k].tf, resp_nt[k], resp_nt[k])
    end
    return resp_nt
end

function forward_helper(
        ::Type{T}, m0, vars, response_trans_utils, params) where {T <: tune_rp_modelType}
    m = tune_rp_modelType(m0.fn_list, T.parameters[2])
    model = from_nt(m, m0)
    resp_nt = to_resp_nt(forward(model, vars, params))
    for k in propertynames(resp_nt)
        broadcast!(response_trans_utils[k].tf, resp_nt[k], resp_nt[k])
    end
    return resp_nt
end

# Only occam uses the following

function zero_abstract(m::mtresponse) where {mtresponse <: MTResponse{
        <:AbstractVector{<:Any}, <:AbstractVector{<:Any}}}
    MTResponse{AbstractVector, AbstractVector}(zero(m.ρₐ), zero(m.ϕ))
end

function inverse(t::mtresponse; abstract=false) where {mtresponse <: MTResponse}
    if abstract
        return MTModel{AbstractArray{<:Any, length(size(t.ρₐ))}, # length(size(...)) => dimensionality
            AbstractArray{<:Any, length(size(t.ρₐ))}}
    else
        vec_type = typeof(t.ρₐ)
        return MTModel{vec_type, vec_type}
    end
end

function inverse(t::rpresponse; abstract=false) where {rpresponse <: RockphyCond}
    if abstract
        return mixing_models{AbstractArray{<:Any, 1}, AbstractArray{<:Any, 1}}
    else
        vec_type = Vector{eltype(t.σ)}
        return mixing_models{vec_type, vec_type}
        # return MTModel{vec_type, vec_type}
    end
end

function forward(t::mtmodel; abstract=false) where {mtmodel <: MTModel} # {<:AbstractVector{<:Any}, <:AbstractVector{<:Any}}}
    if abstract
        return MTResponse{AbstractArray{<:Any, length(size(m.m))},
            AbstractArray{<:Any, length(size(m.m))}}
    else
        vec_type = typeof(t.m)
        return MTResponse{vec_type, vec_type}
    end
end
