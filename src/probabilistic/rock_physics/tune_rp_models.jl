mutable struct tune_rp_modelDistributionType{K, M}
    fn_list::K
    model::Type{M}
end

function tune_rp_modelDistributionType(fn_list, model)
    tune_rp_modelDistributionType(fn_list, typeof(model))
end

mutable struct tune_rp_modelDistribution{K, M, names} <: AbstractRockphyModelDistribution
    fn_list::K
    model::Type{M}
    ps_nt::NamedTuple{names}
end

function (m::tune_rp_modelDistributionType{K, M})(ps_nt) where {K, M}
    return tune_rp_modelDistribution(m.fn_list, m.model, ps_nt)
end

function from_nt(m::T, ps_nt::NamedTuple) where {T <: tune_rp_modelType}
    for fn in m.fn_list
        nt = fn(ps_nt)
        ps_nt = (; ps_nt..., nt...)
    end
    return from_nt(m.model, ps_nt)
end
