@testitem "elastic tests" tags=[:elastic] begin
    using JET
    T = collect(1273.0f0:30:1573.0f0)
    ρ = collect(3300.0f0:100.0f0:4300.0f0)
    ϕ = collect(1.0f-2:1.0f-3:2.0f-2)
    P = 2 .+ zero(T)

    inps = (anharmonic=[T, P, ρ], anharmonic_poro=[T, P, ρ, ϕ], SLB2005=[T, P])

    methods_list = [anharmonic, anharmonic_poro, SLB2005]
    methodsD_list = [
        anharmonicDistribution, anharmonic_poroDistribution, SLB2005Distribution]

    outs = (
        anharmonic=RockphyElastic(
            [7.1367f10, 7.0959f10, 7.0551f10, 7.0143f10, 6.9735f10, 6.9327f10,
                6.8919f10, 6.8511f10, 6.8103f10, 6.7695f10, 6.7287f10],
            [1.1895f11, 1.1827f11, 1.1759f11, 1.1691f11, 1.1623f11, 1.1555f11,
                1.1487f11, 1.1419f11, 1.1351f11, 1.1283f11, 1.1215f11],
            [8.0548f03, 7.9127f03, 7.7764f03, 7.6454f03, 7.5194f03, 7.3981f03,
                7.2811f03, 7.1682f03, 7.0591f03, 6.9537f03, 6.8516f03],
            [4.6504f03, 4.5684f03, 4.4897f03, 4.4141f03, 4.3413f03, 4.2713f03,
                4.2038f03, 4.1386f03, 4.0756f03, 4.0147f03, 3.9558f03]),
        anharmonic_poro=RockphyElastic(
            [6.9053f10, 6.8464f10, 6.7877f10, 6.7293f10, 6.6711f10, 6.6132f10,
                6.5556f10, 6.4982f10, 6.4410f10, 6.3840f10, 6.3272f10],
            [1.1682f11, 1.1596f11, 1.1510f11, 1.1425f11, 1.1341f11, 1.1257f11,
                1.1173f11, 1.1090f11, 1.1007f11, 1.0925f11, 1.0843f11],
            [7.9531f03, 7.8040f03, 7.6610f03, 7.5236f03, 7.3914f03, 7.2641f03,
                7.1414f03, 7.0230f03, 6.9086f03, 6.7980f03, 6.6910f03],
            [4.5744f03, 4.4874f03, 4.4038f03, 4.3235f03, 4.2462f03, 4.1717f03,
                4.0999f03, 4.0306f03, 3.9635f03, 3.8987f03, 3.8359f03]),
        SLB2005=RockphyElastic(zeros(Float32, 11),
            zeros(Float32, 11),
            zeros(Float32, 11),
            [4.4782f03, 4.4669f03, 4.4555f03, 4.4442f03, 4.4328f03, 4.4215f03,
                4.4102f03, 4.3988f03, 4.3875f03, 4.3761f03, 4.3648f03]))

    for i in eachindex(methods_list)
        m = keys(inps)[i]
        model = methods_list[i](inps[m]...)
        out_ = forward(model, [])
        @inferred forward(model, [])
        @test_opt forward(model, [])
        @test_call forward(model, [])
        for k in fieldnames(RockphyElastic)
            @test all(isapprox.(getfield(out_, k), getfield(outs[m], k), rtol=1e-2))
        end
        modelD = methodsD_list[i](inps[m]...)
        @test sample_type(modelD) <: methods_list[i]
        params_ = default_params(typeof(model))
        @inferred forward(model, [], params_)
        @test_opt forward(model, [], params_)
        @test_call forward(model, [], params_)

        tf_ = (; G=no_tf, K=no_tf, Vp=no_tf, Vs=no_tf)
        @inferred forward_helper(sample_type(modelD), to_nt(model), [], tf_, params_)
        @test_opt forward_helper(sample_type(modelD), to_nt(model), [], tf_, params_)
        @test_call forward_helper(sample_type(modelD), to_nt(model), [], tf_, params_)
    end
end
