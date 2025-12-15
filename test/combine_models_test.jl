@testitem "combine models" tags = [:combine_models] begin
    using JET
    m = multi_rp_modelType(SEO3, anharmonic, Nothing, Nothing)
    ps_nt = (; T=[800.0f0, 1000.0f0] .+ 273, P=3.0f0, ρ=3300.0f0, Ch2o_m=1000.0f0, ϕ=0.1f0)
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
end
