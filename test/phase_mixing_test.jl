@testitem "two_phase" tags = [:phase_mixing] begin
    using JET
    m1 = two_phase_modelType(SEO3, Ni2011, HS1962_plus)
    ps_nt = (;
        T = [1200.0f0, 1400.0f0] .+ 273,
        P = 3.0f0,
        ρ = 3300.0f0,
        Ch2o_m = 100.0f0,
        ϕ = 0.1f0,
    )
    model = m1(ps_nt)
    @inferred m1(ps_nt)
    @inferred forward(model, [])

    resp = forward(model, [])
    @test :σ ∈ propertynames(resp)
end

@testitem "multi_phase" tags = [:phase_mixing] begin
    using JET
    m1 = multi_phase_modelType(SEO3, Sifre2014, Zhang2012, HS_minus_multi_phase)
    ps_nt = (;
        T = [1200.0f0] .+ 273,
        Ch2o_ol = 1.0f0,
        Ch2o_m = 1000.0f0,
        Cco2_m = 1000.0f0,
        Ch2o_opx = 100.0f0,
        ϕ = [0.1f0, 0.2f0],
    )
    model = m1(ps_nt)
    @inferred m1(ps_nt)
    @inferred forward(model, [])

    resp = forward(model, [])
    @test :σ ∈ propertynames(resp)
end
