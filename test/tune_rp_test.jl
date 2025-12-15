@testitem "tune rp models" tags=[:tune_rp] begin
    using JET
    T = (800:200:1200) .+ 273.0
    P = [2.0]
    T_solidus = [1000.0 + 273]
    ps_nt = (; T=T, P=P, ϕ=0.1)
    m_type = two_phase_modelType(SEO3, Gaillard2008, HS1962_plus)
    fn_list = [solidus_Hirschmann2000]
    m = tune_rp_modelType(fn_list, m_type)
    model = m(ps_nt)
    @inferred m(ps_nt)
    @report_opt m(ps_nt)
    @report_call m(ps_nt)

    m_dist_ = two_phase_modelDistributionType(SEO3Distribution, Gaillard2008Distribution, HS1962_plus)
    m_dist = tune_rp_modelDistributionType(fn_list, m_dist_)
    modelD = m_dist(ps_nt)

    sample_type(modelD) <: two_phase_modelType
    params_ = (; m1=default_params(SEO3), m2=default_params(Gaillard2008))

    @inferred forward_helper(
        sample_type(modelD), (; ps_nt..., fn_list), [], (; σ=no_tf), params_)
    @test_opt forward_helper(
        sample_type(modelD), (; ps_nt..., fn_list), [], (; σ=no_tf), params_)
    @test_call forward_helper(
        sample_type(modelD), (; ps_nt..., fn_list), [], (; σ=no_tf), params_)
end
