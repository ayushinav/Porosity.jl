@testitem "conductivity tests" tags=[:rp] begin
    using JET
    methods_list = [
        SEO3,
        UHO2014,
        Jones2012,
        Yoshino2009,
        Wang2006,
        Poe2010,

        # melt
        Ni2011,
        Sifre2014,
        Gaillard2008
    ]

    T = collect(1273.0f0:30:1573.0f0)
    Ch2o_ol = collect(2.0f4:2.0f3:4.0f4)
    Ch2o_m = collect(2.0f4:2.0f3:4.0f4)
    Cco2_m = collect(1.0f4:2.0f3:3.0f4)

    inps = (;
        SEO3=[T],
        UHO2014=[T, Ch2o_ol],
        Jones2012=[T, Ch2o_ol],
        Yoshino2009=[T, Ch2o_ol],
        Wang2006=[T, Ch2o_ol],
        Poe2010=[T, Ch2o_ol],

        # melt
        Ni2011=[T, Ch2o_m],
        Sifre2014=[T, Ch2o_m, Cco2_m],
        Gaillard2008=[T]
    )

    methods_list = [
        SEO3,
        UHO2014,
        Jones2012,
        Yoshino2009,
        Wang2006,
        Poe2010,

        # melt
        Ni2011,
        Sifre2014,
        Gaillard2008
    ]

    outs = (;
        SEO3=RockphyCond(
            Float32.([
            -4.0572,
            -3.8956,
            -3.7387,
            -3.5857,
            -3.4355,
            -3.2874,
            -3.1402,
            -2.9933,
            -2.8457,
            -2.6968,
            -2.5463
        ]),
        ),
        UHO2014=RockphyCond(
            Float32.([
            1.2728,
            1.4152,
            1.5459,
            1.6664,
            1.7780,
            1.8818,
            1.9786,
            2.0692,
            2.1542,
            2.2341,
            2.3095
        ]),
        ),
        Jones2012=RockphyCond(
            Float32.([
            0.1550,
            0.2774,
            0.3920,
            0.4996,
            0.6011,
            0.6971,
            0.7880,
            0.8744,
            0.9567,
            1.0351,
            1.1099
        ]),
        ),
        Yoshino2009=RockphyCond(
            Float32.([
            -0.6429,
            -0.5108,
            -0.3878,
            -0.2729,
            -0.1650,
            -0.0634,
            0.0325,
            0.1232,
            0.2093,
            0.2912,
            0.3691
        ]),
        ),
        Wang2006=RockphyCond(
            Float32.([
            -0.3832,
            -0.2753,
            -0.1734,
            -0.0768,
            0.0150,
            0.1023,
            0.1857,
            0.2653,
            0.3415,
            0.4144,
            0.4844
        ]),
        ),
        Poe2010=RockphyCond(
            Float32.([
            3.4473,
            3.6441,
            3.8203,
            3.9789,
            4.1224,
            4.2527,
            4.3714,
            4.4800,
            4.5796,
            4.6711,
            4.7554
        ]),
        ),
        Ni2011=RockphyCond(
            Float32.([
            -2.3579,
            -1.3975,
            -0.7500,
            -0.2847,
            0.0652,
            0.3375,
            0.5552,
            0.7329,
            0.8807,
            1.0053,
            1.1117
        ]),
        ),
        Sifre2014=RockphyCond(
            Float32.([
            -0.0128,
            0.1570,
            0.3119,
            0.4537,
            0.5837,
            0.7031,
            0.8131,
            0.9145,
            1.0084,
            1.0954,
            1.1763
        ]),
        ),
        Gaillard2008=RockphyCond(
            Float32.([
            2.2276,
            2.2577,
            2.2865,
            2.3140,
            2.3403,
            2.3655,
            2.3897,
            2.4129,
            2.4352,
            2.4566,
            2.4772
        ]),
        )
    )

    for i in eachindex(methods_list)
        m = keys(inps)[i]
        model = methods_list[i](inps[m]...)
        out_ = forward(model, [])
        @inferred forward(model, [])
        @report_call forward(model, [])
        for k in fieldnames(RockphyCond)
            @test all(isapprox.(getfield(out_, k), getfield(outs[m], k), rtol=1e-2))
        end
    end
end

@testitem "elastic tests" tags=[:rp] begin
    using JET
    T = collect(1273.0f0:30:1573.0f0)
    ρ = collect(3300.0f0:100.0f0:4300.0f0)
    ϕ = collect(1.0f-2:1.0f-3:2.0f-2)
    P = 2 .+ zero(T)

    inps = (anharmonic=[T, P, ρ], anharmonic_poro=[T, P, ρ, ϕ], SLB2005=[T, P])

    outs = (
        anharmonic=MT.RockphyElastic(
            [
                7.1367f10,
                7.0959f10,
                7.0551f10,
                7.0143f10,
                6.9735f10,
                6.9327f10,
                6.8919f10,
                6.8511f10,
                6.8103f10,
                6.7695f10,
                6.7287f10
            ],
            [
                1.1895f11,
                1.1827f11,
                1.1759f11,
                1.1691f11,
                1.1623f11,
                1.1555f11,
                1.1487f11,
                1.1419f11,
                1.1351f11,
                1.1283f11,
                1.1215f11
            ],
            [
                8.0548f03,
                7.9127f03,
                7.7764f03,
                7.6454f03,
                7.5194f03,
                7.3981f03,
                7.2811f03,
                7.1682f03,
                7.0591f03,
                6.9537f03,
                6.8516f03
            ],
            [
                4.6504f03,
                4.5684f03,
                4.4897f03,
                4.4141f03,
                4.3413f03,
                4.2713f03,
                4.2038f03,
                4.1386f03,
                4.0756f03,
                4.0147f03,
                3.9558f03
            ]
        ),
        anharmonic_poro=RockphyElastic(
            [
                6.9053f10,
                6.8464f10,
                6.7877f10,
                6.7293f10,
                6.6711f10,
                6.6132f10,
                6.5556f10,
                6.4982f10,
                6.4410f10,
                6.3840f10,
                6.3272f10
            ],
            [
                1.1682f11,
                1.1596f11,
                1.1510f11,
                1.1425f11,
                1.1341f11,
                1.1257f11,
                1.1173f11,
                1.1090f11,
                1.1007f11,
                1.0925f11,
                1.0843f11
            ],
            [
                7.9531f03,
                7.8040f03,
                7.6610f03,
                7.5236f03,
                7.3914f03,
                7.2641f03,
                7.1414f03,
                7.0230f03,
                6.9086f03,
                6.7980f03,
                6.6910f03
            ],
            [
                4.5744f03,
                4.4874f03,
                4.4038f03,
                4.3235f03,
                4.2462f03,
                4.1717f03,
                4.0999f03,
                4.0306f03,
                3.9635f03,
                3.8987f03,
                3.8359f03
            ]
        ),
        SLB2005=RockphyElastic(
            zeros(Float32, 11),
            zeros(Float32, 11),
            zeros(Float32, 11),
            [
                4.4782f03,
                4.4669f03,
                4.4555f03,
                4.4442f03,
                4.4328f03,
                4.4215f03,
                4.4102f03,
                4.3988f03,
                4.3875f03,
                4.3761f03,
                4.3648f03
            ]
        )
    )

    methods_list = [anharmonic, anharmonic_poro, SLB2005]

    for i in eachindex(methods_list)
        m = keys(inps)[i]
        model = methods_list[i](inps[m]...)
        out_ = forward(model, [])
        @inferred forward(model, [])
        @report_call forward(model, [])
        for k in fieldnames(RockphyElastic)
            @test all(isapprox.(getfield(out_, k), getfield(outs[m], k), rtol=1e-2))
        end
    end
end

@testitem "viscosity tests" tags=[:rp] begin
    using JET
    T = collect(1073.0f0:30:1373.0f0)
    P = 2 .+ zero(T)
    dg = collect(3.0f0:4.0f-1:7.0f0)
    σ = collect(7.5f0:0.5f0:12.5f0) .* 1.0f-3
    ϕ = collect(1.0f-2:1.0f-3:2.0f-2)
    T_solidus = 1473 .+ zero(T)

    inps = (
        HZK2011=[T, P, dg, σ, ϕ],
        HK2003=[T, P, dg, σ, ϕ, zero(ϕ)],
        xfit_premelt=[T, P, dg, σ, ϕ, T_solidus]
    )

    methods_list = [HZK2011, HK2003, xfit_premelt]

    outs = (
        HZK2011=RockphyViscous(
            Float32[
                1.0f-12,
                2.0f-12,
                5.0f-12,
                1.2f-11,
                2.8f-11,
                6.2f-11,
                1.36f-10,
                2.86f-10,
                5.86f-10,
                1.1710001f-9,
                2.283f-9
            ],
            [
                8.9619f18,
                3.8154f18,
                1.6603f18,
                7.410f17,
                3.397f17,
                1.601f17,
                7.75f16,
                3.85f16,
                1.96f16,
                1.03f16,
                5.5f15
            ]
        ),
        HK2003=RockphyViscous(
            Float32[
                3.0f-11,
                7.999999f-11,
                1.8999999f-10,
                4.5999998f-10,
                1.05f-9,
                2.3499998f-9,
                5.09f-9,
                1.0709999f-8,
                2.191f-8,
                4.3609997f-8,
                8.4579995f-8
            ],
            [
                2.3787f17,
                1.0128f17,
                4.408f16,
                1.968f16,
                9.03f15,
                4.26f15,
                2.06f15,
                1.03f15,
                5.2f14,
                2.8f14,
                1.5f14
            ]
        ),
        xfit_premelt=RockphyViscous(
            zero(T),
            [
                7.6341f18,
                2.5851f18,
                9.069f17,
                3.305f17,
                1.251f17,
                4.92f16,
                2.01f16,
                8.5f15,
                3.7f15,
                1.7f15,
                8.0f14
            ]
        )
    )

    for i in eachindex(methods_list)
        m = keys(inps)[i]
        model = methods_list[i](inps[m]...)
        out_ = forward(model, [])
        @inferred forward(model, [])
        @report_call forward(model, [])
        (out_.ϵ_rate .- outs[m].ϵ_rate) ./ outs[m].ϵ_rate
        (out_.η .- outs[m].η) ./ outs[m].η
        for k in fieldnames(RockphyViscous)
            @test all(
                isapprox.(
                log10.(getfield(out_, k)),
                log10.(getfield(outs[m], k)),
                rtol=1.0f-2
            ),
            )
        end
    end
end

@testitem "anelastic tests" tags=[:rp] begin
    using JET
    T = collect(1073.0f0:30:1373.0f0)
    P = 2 .+ zero(T)
    dg = collect(3.0f0:4.0f-1:7.0f0)
    σ = collect(7.5f0:0.5f0:12.5f0) .* 1.0f-3
    ϕ = collect(1.0f-2:1.0f-3:2.0f-2)
    T_solidus = 1473 .+ zero(T)
    f = [1.0f0] #10f0 .^ collect(-10:1:0)
    ρ = collect(3300.0f0:100.0f0:4300.0f0)

    inps = (
        andrade_psp=[T, P, dg, σ, ϕ, ρ, f'],
        eburgers_psp=[T, P, dg, σ, ϕ, ρ, zero(ϕ), T_solidus, f'],
        premelt_anelastic=[T, P, dg, σ, ϕ, ρ, zero(ϕ), T_solidus, f'],
        xfit_mxw=[T, P, dg, σ, ϕ, ρ, zero(ϕ), T_solidus, f'],
        andrade_analytical=[T, P, dg, σ, ϕ, ρ, zero(ϕ), T_solidus, f']
    )

    methods_list = [
        andrade_psp, eburgers_psp, premelt_anelastic, xfit_mxw, andrade_analytical]

    outs = (
        andrade_psp=RockphyAnelastic(
            [
                1.352f-11,
                1.361f-11,
                1.369f-11,
                1.378f-11,
                1.388f-11,
                1.398f-11,
                1.408f-11,
                1.419f-11,
                1.431f-11,
                1.444f-11,
                1.458f-11
            ],
            [
                1.40f-14,
                1.89f-14,
                2.51f-14,
                3.30f-14,
                4.29f-14,
                5.51f-14,
                7.00f-14,
                8.80f-14,
                1.098f-13,
                1.356f-13,
                1.663f-13
            ],
            [
                1.039f-3,
                1.39f-3,
                1.836f-3,
                2.397f-3,
                3.090f-3,
                3.9398f-3,
                4.9691f-3,
                6.20326f-3,
                7.669f-3,
                9.394f-3,
                1.141f-2
            ],
            [
                7.3952f10,
                7.3499f10,
                7.3035f10,
                7.2557f10,
                7.2062f10,
                7.1549f10,
                7.1014f10,
                7.0455f10,
                6.9869f10,
                6.9253f10,
                6.8604f10
            ],
            [
                4.7339f03,
                4.6495f03,
                4.5681f03,
                4.4894f03,
                4.4132f03,
                4.3392f03,
                4.2672f03,
                4.1969f03,
                4.1281f03,
                4.0606f03,
                3.9943f03
            ],
            [
                4.7339f03,
                4.6495f03,
                4.5681f03,
                4.4894f03,
                4.4132f03,
                4.3392f03,
                4.2672f03,
                4.1969f03,
                4.1281f03,
                4.0606f03,
                3.9943f03
            ]
        ),
        eburgers_psp=RockphyAnelastic(
            [
                1.353f-11,
                1.362f-11,
                1.372f-11,
                1.382f-11,
                1.393f-11,
                1.405f-11,
                1.417f-11,
                1.430f-11,
                1.443f-11,
                1.458f-11,
                1.473f-11
            ],
            [
                3.06f-14,
                3.92f-14,
                4.90f-14,
                6.05f-14,
                7.40f-14,
                8.97f-14,
                1.078f-13,
                1.287f-13,
                1.525f-13,
                1.795f-13,
                2.099f-13
            ],
            [
                2.3f-3,
                2.9f-3,
                3.6f-3,
                4.4f-3,
                5.3f-3,
                6.4f-3,
                7.6f-3,
                9.0f-3,
                1.06f-2,
                1.23f-2,
                1.42f-2
            ],
            [
                7.3907f10,
                7.3401f10,
                7.2878f10,
                7.2335f10,
                7.1771f10,
                7.1185f10,
                7.0576f10,
                6.9941f10,
                6.9279f10,
                6.8590f10,
                6.7872f10
            ],
            [
                4.7325f03,
                4.6464f03,
                4.5631f03,
                4.4825f03,
                4.4043f03,
                4.3282f03,
                4.2540f03,
                4.1815f03,
                4.1106f03,
                4.0412f03,
                3.9729f03
            ],
            [
                4.7325f03,
                4.6464f03,
                4.5631f03,
                4.4825f03,
                4.4043f03,
                4.3282f03,
                4.2540f03,
                4.1815f03,
                4.1106f03,
                4.0412f03,
                3.9729f03
            ]
        ),
        premelt_anelastic=RockphyAnelastic(
            [
                1.351f-11,
                1.360f-11,
                1.370f-11,
                1.380f-11,
                1.393f-11,
                1.407f-11,
                1.423f-11,
                1.442f-11,
                1.464f-11,
                1.491f-11,
                1.576f-11
            ],
            [
                1.28f-14,
                2.23f-14,
                3.73f-14,
                5.92f-14,
                8.89f-14,
                1.269f-13,
                1.728f-13,
                2.259f-13,
                2.854f-13,
                3.646f-13,
                6.147f-13
            ],
            [
                9.4591f-4,
                1.6425f-3,
                2.7f-3,
                4.3f-3,
                6.4f-3,
                9.0f-3,
                1.21f-2,
                1.57f-2,
                1.95f-2,
                2.44f-2,
                3.90f-2
            ],
            [
                7.3998f10,
                7.3527f10,
                7.3016f10,
                7.2448f10,
                7.1808f10,
                7.1081f10,
                7.0256f10,
                6.9328f10,
                6.8299f10,
                6.7042f10,
                6.3399f10
            ],
            [
                4.7354f03,
                4.6503f03,
                4.5675f03,
                4.4861f03,
                4.4055f03,
                4.3251f03,
                4.2445f03,
                4.1634f03,
                4.0818f03,
                3.9959f03,
                3.8412f03
            ],
            [
                4.7354f03,
                4.6503f03,
                4.5675f03,
                4.4861f03,
                4.4055f03,
                4.3251f03,
                4.2445f03,
                4.1634f03,
                4.0818f03,
                3.9959f03,
                3.8412f03
            ]
        ),
        xfit_mxw=RockphyAnelastic(
            [
                1.412f-11,
                1.432f-11,
                1.453f-11,
                1.475f-11,
                1.498f-11,
                1.522f-11,
                1.548f-11,
                1.574f-11,
                1.603f-11,
                1.633f-11,
                1.665f-11
            ],
            [
                1.713f-13,
                1.871f-13,
                2.051f-13,
                2.260f-13,
                2.501f-13,
                2.779f-13,
                3.102f-13,
                3.476f-13,
                3.909f-13,
                4.410f-13,
                4.991f-13
            ],
            [
                0.0121,
                0.0131,
                0.0141,
                0.0153,
                0.0167,
                0.0183,
                0.0200,
                0.0221,
                0.0244,
                0.0270,
                0.0300
            ],
            [
                7.0802f10,
                6.9807f10,
                6.8797f10,
                6.7773f10,
                6.6733f10,
                6.5676f10,
                6.4598f10,
                6.3498f10,
                6.2372f10,
                6.1218f10,
                6.0033f10
            ],
            [
                4.6322f03,
                4.5314f03,
                4.4338f03,
                4.3391f03,
                4.2472f03,
                4.1576f03,
                4.0702f03,
                3.9848f03,
                3.9009f03,
                3.8185f03,
                3.7373f03
            ],
            [
                4.6322f03,
                4.5314f03,
                4.4338f03,
                4.3391f03,
                4.2472f03,
                4.1576f03,
                4.0702f03,
                3.9848f03,
                3.9009f03,
                3.8185f03,
                3.7373f03
            ]
        ),
        andrade_analytical=RockphyAnelastic(
            [
                1.3500f-11,
                1.3570f-11,
                1.3650f-11,
                1.3720f-11,
                1.3800f-11,
                1.3880f-11,
                1.3960f-11,
                1.4040f-11,
                1.4120f-11,
                1.4200f-11,
                1.4280f-11
            ],
            [
                3.27f-16,
                3.2999999f-16,
                3.34f-16,
                3.4f-16,
                3.52f-16,
                3.7299998f-16,
                4.1499998f-16,
                4.95f-16,
                6.45f-16,
                9.22f-16,
                1.422f-15
            ],
            [
                2.4240f-5,
                2.4310f-5,
                2.4460f-5,
                2.4780f-5,
                2.5470f-5,
                2.6890f-5,
                2.9720f-5,
                3.5230f-5,
                4.5650f-5,
                6.4890f-5,
                9.9510f-5
            ],
            [
                7.4084f10,
                7.3676f10,
                7.3268f10,
                7.2860f10,
                7.2452f10,
                7.2044f10,
                7.1636f10,
                7.1228f10,
                7.0820f10,
                7.0412f10,
                7.0004f10
            ],
            [
                4.7381f03,
                4.6550f03,
                4.5753f03,
                4.4988f03,
                4.4251f03,
                4.3542f03,
                4.2858f03,
                4.2198f03,
                4.1561f03,
                4.0945f03,
                4.0349f03
            ],
            [
                4.7381f03,
                4.6550f03,
                4.5753f03,
                4.4988f03,
                4.4251f03,
                4.3542f03,
                4.2858f03,
                4.2198f03,
                4.1561f03,
                4.0945f03,
                4.0349f03
            ]
        )
    )

    for i in eachindex(methods_list)
        m = keys(inps)[i]
        model = methods_list[i](inps[m]...)
        out_ = forward(model, [])
        # @inferred forward(model, []) (cannot infer using quadgk)
        @report_call forward(model, [])
        for k in fieldnames(RockphyAnelastic)
            @test all(
                isapprox.(
                log10.(getfield(out_, k)),
                log10.(getfield(outs[m], k)),
                rtol=1.5f-2
            ),
            )
        end
    end
end

@testitem "two_phase" tags=[:rp] begin
    using JET
    m1 = two_phase_modelType(SEO3, Ni2011, HS1962_plus())
    ps_nt = (; T=[1200.0f0, 1400.0f0] .+ 273, P=3.0f0, ρ=3300.0f0, Ch2o_m=100.0f0, ϕ=0.1f0)
    model = m1(ps_nt)
    @inferred m1(ps_nt)
    @inferred forward(model, [])
end

@testitem "multi_phase" tags=[:rp] begin
    using JET
    m1 = multi_phase_modelType(SEO3, Sifre2014, Zhang2012, HS_minus_multi_phase())
    ps_nt = (; T=[1200.0f0] .+ 273, Ch2o_ol=1.0f0, Ch2o_m=1000.0f0,
        Cco2_m=1000.0f0, Ch2o_opx=100.0f0, ϕ=[0.1f0, 0.2f0])
    model = m1(ps_nt)
    @inferred m1(ps_nt)
    @inferred forward(model, [])
end

@testitem "combine models" tags=[:rp] begin
    using JET
    m = multi_rp_modelType(SEO3, anharmonic, Nothing, Nothing)
    ps_nt = (; T=[800.0f0, 1000.0f0] .+ 273, P=3.0f0, ρ=3300.0f0, Ch2o_m=1000.0f0, ϕ=0.1f0)
    model = m(ps_nt)
    resp = MT.to_resp_nt(forward(model, []))
    @inferred forward(model, [])

    resp_SEO3 = forward(model.cond, [])
    resp_anharmonic = forward(model.elastic, [])

    resp_check = merge(MT.to_nt.([resp_SEO3, resp_anharmonic])...)

    @test resp == resp_check
end
