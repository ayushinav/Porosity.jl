# NamedTuple manipulation

function to_resp_nt(d::T) where {T <: multi_rp_response}
    return merge(to_nt(d.cond), to_nt(d.elastic), to_nt(d.visc), to_nt(d.anelastic))
end

# forward manipulation

function SubsurfaceCore.forward_helper(
        m::Type{T}, m0, vars, response_trans_utils, params) where {T <: two_phase_modelType}
    model = from_nt(m, m0)
    resp_nt = to_resp_nt(forward(model, vars, params))
    for k in propertynames(resp_nt)
        broadcast!(response_trans_utils[k].tf, resp_nt[k], resp_nt[k])
    end
    return resp_nt
end

function SubsurfaceCore.forward_helper(m::Type{T}, m0, vars, response_trans_utils,
        params) where {T <: multi_phase_modelType}
    model = from_nt(m, m0)
    resp_nt = to_resp_nt(forward(model, vars, params))
    for k in propertynames(resp_nt)
        broadcast!(response_trans_utils[k].tf, resp_nt[k], resp_nt[k])
    end
    return resp_nt
end

function SubsurfaceCore.forward_helper(
        m::Type{T}, m0, vars, response_trans_utils, params) where {T <: multi_rp_modelType}
    model = from_nt(m, m0)
    resp_nt = to_resp_nt(forward(model, vars, params))
    for k in propertynames(resp_nt)
        broadcast!(response_trans_utils[k].tf, resp_nt[k], resp_nt[k])
    end
    return resp_nt
end

function SubsurfaceCore.forward_helper(
        ::Type{T}, m0, vars, response_trans_utils, params) where {T <: tune_rp_modelType}
    m = tune_rp_modelType(m0.fn_list, T.parameters[2])
    model = from_nt(m, m0)
    resp_nt = to_resp_nt(forward(model, vars, params))
    for k in propertynames(resp_nt)
        broadcast!(response_trans_utils[k].tf, resp_nt[k], resp_nt[k])
    end
    return resp_nt
end
