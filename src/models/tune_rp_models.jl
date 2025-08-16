mutable struct tune_rp_modelType{K, M} <: AbstractRockphyModel
    fn_list::K
    model::Type{M}
end

tune_rp_modelType(fn_list, model) = tune_rp_modelType(fn_list, typeof(model))

function (m::tune_rp_modelType{K, M})(ps_nt) where {K, M}
    for fn in m.fn_list
        nt = fn(ps_nt)
        ps_nt = (; ps_nt..., nt...)
    end
    return from_nt(M, ps_nt)
end

function (m::tune_rp_modelType{K, M})(ps_nt) where {K, M <: two_phase_modelType}
    for fn in m.fn_list
        nt = fn(ps_nt)
        ps_nt = (; ps_nt..., nt...)
    end
    return from_nt(m.model, ps_nt)
end

default_params(::Type{tune_rp_modelType{K, M}}) where {K, M} = default_params(M)
