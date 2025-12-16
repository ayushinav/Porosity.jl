@testitem "two_phase" tags=[:phase_mixing] begin
    using JET
    m1 = two_phase_modelType(SEO3, Gaillard2008, HS1962_plus)
    ps_nt = (; T=[1200.0f0, 1400.0f0] .+ 273, P=3.0f0, ρ=3300.0f0, Ch2o_m=100.0f0, ϕ=0.1f0)
    model = m1(ps_nt)
    @inferred m1(ps_nt)
    @test_opt m1(ps_nt)
    @test_call m1(ps_nt)

    @inferred forward(model, [])
    @test_opt forward(model, [])
    @test_call forward(model, [])
    resp = forward(model, [])
    @test :σ ∈ propertynames(resp)

    m1dist = two_phase_modelDistributionType(
        SEO3Distribution, Gaillard2008Distribution, HS1962_plus)
    modelD = m1dist(ps_nt)
    @test sample_type(modelD) <: two_phase_modelType

    params_ = (; m1=default_params(SEO3), m2=default_params(Gaillard2008))

    @inferred forward_helper(sample_type(modelD), ps_nt, [], (; σ=no_tf), params_)
    @test_opt forward_helper(sample_type(modelD), ps_nt, [], (; σ=no_tf), params_)
    @test_call forward_helper(sample_type(modelD), ps_nt, [], (; σ=no_tf), params_)
end

@testitem "multi_phase" tags=[:phase_mixing] begin
    using JET
    m1 = multi_phase_modelType(SEO3, Sifre2014, Zhang2012, HS_minus_multi_phase)
    ps_nt = (; T=[1200.0f0] .+ 273, Ch2o_ol=[1.0f0], Ch2o_m=[1000.0f0],
        Cco2_m=[1000.0f0], Ch2o_opx=[100.0f0], ϕ=[0.1f0, 0.2f0])
    model = m1(ps_nt)
    @inferred m1(ps_nt)
    @test_opt m1(ps_nt)
    @test_call m1(ps_nt)

    @inferred forward(model, [])
    @test_opt forward(model, [])
    @test_call forward(model, [])
    resp = forward(model, [])

    @test :σ ∈ propertynames(resp)

    m1dist = multi_phase_modelDistributionType(SEO3Distribution, Sifre2014Distribution,
        Zhang2012Distribution, HS_minus_multi_phase)
    modelD = m1dist(ps_nt)
    @test sample_type(modelD) <: multi_phase_modelType

    params_ = (; m1=default_params(SEO3), m2=default_params(Sifre2014),
        m3=default_params(Zhang2012))

    @inferred forward_helper(sample_type(modelD), ps_nt, [], (; σ=no_tf), params_)
    @test_opt forward_helper(sample_type(modelD), ps_nt, [], (; σ=no_tf), params_)
end
