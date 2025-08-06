for m in subtypes(AbstractGeophyModel)
    mstring = string(m)
    dstring = mstring * "Distribution"
    eval(Meta.parse("SubsurfaceCore.sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("SubsurfaceCore.sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end

for m in subtypes(AbstractCondModel)
    mstring = string(m)
    dstring = mstring * "Distribution"
    eval(Meta.parse("SubsurfaceCore.sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("SubsurfaceCore.sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end

for m in subtypes(AbstractElasticModel)
    mstring = string(m)
    dstring = mstring * "Distribution"
    eval(Meta.parse("SubsurfaceCore.sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("SubsurfaceCore.sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end

for m in subtypes(AbstractViscousModel)
    mstring = string(m)
    dstring = mstring * "Distribution"
    eval(Meta.parse("SubsurfaceCore.sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("SubsurfaceCore.sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end

for m in subtypes(AbstractAnelasticModel)
    mstring = string(m)
    dstring = mstring * "Distribution"
    eval(Meta.parse("SubsurfaceCore.sample_type(d::$dstring) = $mstring"))
    eval(Meta.parse("SubsurfaceCore.sample_type(::Type{T}) where {T <:$dstring} = $mstring"))
end

function SubsurfaceCore.sample_type(d::two_phase_modelDistribution{
        V, T1, T2, M}) where {V, T1, T2, M}
    t1 = sample_type(T1)
    t2 = sample_type(T2)
    m = M
    two_phase_modelType{t1, t2, m}
end

function SubsurfaceCore.sample_type(d::multi_phase_modelDistribution{
        V, T1, T2, T3, T4, T5, T6, T7, T8, M}) where {V, T1, T2, T3, T4, T5, T6, T7, T8, M}
    t1 = sample_type(T1)
    t2 = sample_type(T2)
    t3 = sample_type(T3)
    t4 = sample_type(T4)
    t5 = sample_type(T5)
    t6 = sample_type(T6)
    t7 = sample_type(T7)
    t8 = sample_type(T8)
    m = M
    multi_phase_modelType{t1, t2, t3, t4, t5, t6, t7, t8, m}
end

function SubsurfaceCore.sample_type(::Type{two_phase_modelDistribution{
        V, T1, T2, M}}) where {V, T1, T2, M}
    # needed for higher models
    t1 = sample_type(T1)
    t2 = sample_type(T2)
    m = sample_type(M)
    two_phase_modelType{t1, t2, m}
end

function SubsurfaceCore.sample_type(::Type{multi_phase_modelDistribution{
        V, T1, T2, T3, T4, T5, T6, T7, T8, M}}) where {V, T1, T2, T3, T4, T5, T6, T7, T8, M}
    # needed for higher models
    t1 = sample_type(T1)
    t2 = sample_type(T2)
    t3 = sample_type(T3)
    t4 = sample_type(T4)
    t5 = sample_type(T5)
    t6 = sample_type(T6)
    t7 = sample_type(T7)
    t8 = sample_type(T8)
    m = sample_type(M)
    multi_phase_modelType{t1, t2, t3, t4, t5, t6, t7, t8, m}
end

function SubsurfaceCore.sample_type(::multi_rp_modelDistribution{
        T1, T2, T3, T4}) where {T1, T2, T3, T4}
    multi_rp_modelType{sample_type(T1), sample_type(T2), sample_type(T3), sample_type(T4)}
end

function SubsurfaceCore.sample_type(d::tune_rp_modelDistribution{K, M}) where {K, M}
    m_ = sample_type(d.model)
    tune_rp_modelType{Vector{Function}, m_}
end

# NamedTuple manipulation

function SubsurfaceCore.to_dist_nt(d::T) where {T <: two_phase_modelDistribution}
    m1 = to_nt(d.m1)
    m2 = to_nt(d.m2)
    mix = to_nt(d.mix)
    return (; ϕ=d.ϕ, m1..., m2..., mix...)
end

function SubsurfaceCore.to_dist_nt(d::T) where {T <: multi_phase_modelDistribution}
    m1 = to_nt(d.m1)
    m2 = to_nt(d.m2)
    m3 = to_nt(d.m3)
    m4 = to_nt(d.m4)
    m5 = to_nt(d.m5)
    m6 = to_nt(d.m6)
    m7 = to_nt(d.m7)
    m8 = to_nt(d.m8)
    mix = to_nt(d.mix)
    return (; ϕ=d.ϕ, m1..., m2..., m3..., m4..., m5..., m6..., m7..., m8..., mix...)
end

function SubsurfaceCore.to_dist_nt(d::T) where {T <: multi_rp_modelDistribution}
    return merge(to_dist_nt(d.cond), to_dist_nt(d.elastic),
        to_dist_nt(d.visc), to_dist_nt(d.anelastic))
end

function SubsurfaceCore.to_dist_nt(d::T) where {T <: multi_rp_responseDistribution}
    return merge(to_nt(d.cond), to_nt(d.elastic), to_nt(d.visc), to_nt(d.anelastic))
end

function SubsurfaceCore.to_dist_nt(d::tune_rp_modelDistribution)
    return (; d.ps_nt..., fn_list=d.fn_list)
end
