@testitem "viscosity tests" tags = [:visc] begin
    using JET
    T = collect(1073.0f0:30:1373.0f0)
    P = 2 .+ zero(T)
    dg = collect(3.0f0:4.0f-1:7.0f0)
    σ = collect(7.5f0:0.5f0:12.5f0) .* 1.0f-3
    ϕ = collect(1.0f-2:1.0f-3:2.0f-2)
    T_solidus = 1473 .+ zero(T)

    inps = (HZK2011=[T, P, dg, σ, ϕ], HK2003=[T, P, dg, σ, ϕ, zero(ϕ)],
        xfit_premelt=[T, P, dg, σ, ϕ, T_solidus])

    methods_list = [HZK2011, HK2003, xfit_premelt]
    methodsD_list = [HZK2011Distribution, HK2003Distribution, xfit_premeltDistribution]

    outs = (
        HZK2011=RockphyViscous(
            Float32[1.0f-12, 2.0f-12, 5.0f-12, 1.2f-11, 2.8f-11, 6.2f-11, 1.36f-10, 2.86f-10, 5.86f-10, 1.1710001f-9, 2.283f-9],
            [8.9619f18, 3.8154f18, 1.6603f18, 7.410f17, 3.397f17,
                1.601f17, 7.75f16, 3.85f16, 1.96f16, 1.03f16, 5.5f15]),
        HK2003=RockphyViscous(
            Float32[3.0f-11, 7.999999f-11, 1.8999999f-10, 4.5999998f-10, 1.05f-9, 2.3499998f-9, 5.09f-9, 1.0709999f-8, 2.191f-8, 4.3609997f-8, 8.4579995f-8],
            [2.3787f17, 1.0128f17, 4.408f16, 1.968f16, 9.03f15,
                4.26f15, 2.06f15, 1.03f15, 5.2f14, 2.8f14, 1.5f14]),
        xfit_premelt=RockphyViscous(zero(T),
            [7.6341f18, 2.5851f18, 9.069f17, 3.305f17, 1.251f17,
                4.92f16, 2.01f16, 8.5f15, 3.7f15, 1.7f15, 8.0f14]))

    for i in eachindex(methods_list)
        m = keys(inps)[i]
        model = methods_list[i](inps[m]...)
        out_ = forward(model, [])
        @inferred forward(model, [])
        @test_opt forward(model, [])
        @test_call forward(model, [])
        for k in fieldnames(RockphyViscous)
            @test all(isapprox.(log10.(getfield(out_, k)), log10.(getfield(outs[m], k)), rtol=1.0f-2),)
        end
        modelD = methodsD_list[i](inps[m]...)
        @test sample_type(modelD) <: methods_list[i]
        params_ = default_params(typeof(model))
        @inferred forward(model, [], params_)
        @test_opt forward(model, [], params_)
        @test_call forward(model, [], params_)

        @inferred forward_helper(
            sample_type(modelD), to_nt(model), [], (; ϵ_rate=no_tf, η=no_tf), params_)
        @test_opt forward_helper(
            sample_type(modelD), to_nt(model), [], (; ϵ_rate=no_tf, η=no_tf), params_)
        @test_call forward_helper(
            sample_type(modelD), to_nt(model), [], (; ϵ_rate=no_tf, η=no_tf), params_)
    end
end
