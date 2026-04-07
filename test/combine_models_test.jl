@testitem "combine models" tags=[:combine_models] begin
    using JET
    m = multi_rp_modelType(SEO3, anharmonic, Nothing, Nothing)
    ps_nt = (;
        T=[800.0f0, 1000.0f0] .+ 273, P=[3.0f0], ρ=[3300.0f0], Ch2o_m=1000.0f0, ϕ=0.1f0)
    model = m(ps_nt)
    @inferred m(ps_nt)
    @test_opt m(ps_nt)
    @test_call m(ps_nt)

    resp = forward(model, [])
    @inferred forward(model, [])
    @test_opt forward(model, [])
    @test_call forward(model, [])

    resp_SEO3 = forward(model.cond, [])
    resp_anharmonic = forward(model.elastic, [])

    @test resp.cond.σ == resp_SEO3.σ
    @test resp.elastic.Vp == resp_anharmonic.Vp
    @test resp.elastic.Vs == resp_anharmonic.Vs

    m_dist = multi_rp_modelDistributionType(
        SEO3Distribution, anharmonicDistribution, Nothing, Nothing)
    modelD = m_dist(ps_nt)

    sample_type(modelD) <: multi_rp_modelType

    params_ = (cond=default_params(SEO3), elastic=default_params(anharmonic),
        visc=default_params(Nothing), anelastic=default_params(Nothing))

    tf_ = (; σ=no_tf, G=no_tf, K=no_tf, Vp=no_tf, Vs=no_tf)

    @inferred forward_helper(sample_type(modelD), ps_nt, [], tf_, params_)
    @test_opt forward_helper(sample_type(modelD), ps_nt, [], tf_, params_)
    @test_call forward_helper(sample_type(modelD), ps_nt, [], tf_, params_)
end
