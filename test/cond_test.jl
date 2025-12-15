@testitem "conductivity tests" tags = [:cond] begin
    using JET

    T = collect(1273.0f0:30:1573.0f0)
    Ch2o_ol = collect(2.0f4:2.0f3:4.0f4)
    Ch2o_m = collect(2.0f4:2.0f3:4.0f4)
    Cco2_m = collect(1.0f4:2.0f3:3.0f4)

    inps = (; SEO3=[T], UHO2014=[T, Ch2o_ol], Jones2012=[T, Ch2o_ol],
        Yoshino2009=[T, Ch2o_ol], Wang2006=[T, Ch2o_ol], Poe2010=[T, Ch2o_ol],

        # melt
        Ni2011=[T, Ch2o_m], Sifre2014=[T, Ch2o_m, Cco2_m], Gaillard2008=[T])

    methods_list = [SEO3, UHO2014, Jones2012, Yoshino2009, Wang2006, Poe2010,

        # melt
        Ni2011, Sifre2014, Gaillard2008]

    methodsD_list = [SEO3Distribution, UHO2014Distribution, Jones2012Distribution,
        Yoshino2009Distribution, Wang2006Distribution, Poe2010Distribution,

        # melt
        Ni2011Distribution, Sifre2014Distribution, Gaillard2008Distribution]

    methods_list = [SEO3, UHO2014, Jones2012, Yoshino2009, Wang2006, Poe2010,

        # melt
        Ni2011, Sifre2014, Gaillard2008]

    outs = (;
        SEO3=RockphyCond(Float32.([-4.0572, -3.8956, -3.7387, -3.5857, -3.4355, -3.2874,
            -3.1402, -2.9933, -2.8457, -2.6968, -2.5463]),),
        UHO2014=RockphyCond(Float32.([1.2728, 1.4152, 1.5459, 1.6664, 1.7780, 1.8818,
            1.9786, 2.0692, 2.1542, 2.2341, 2.3095]),),
        Jones2012=RockphyCond(Float32.([0.1550, 0.2774, 0.3920, 0.4996, 0.6011, 0.6971,
            0.7880, 0.8744, 0.9567, 1.0351, 1.1099]),),
        Yoshino2009=RockphyCond(Float32.([-0.6429, -0.5108, -0.3878, -0.2729, -0.1650,
            -0.0634, 0.0325, 0.1232, 0.2093, 0.2912, 0.3691]),),
        Wang2006=RockphyCond(Float32.([-0.3832, -0.2753, -0.1734, -0.0768, 0.0150,
            0.1023, 0.1857, 0.2653, 0.3415, 0.4144, 0.4844]),),
        Poe2010=RockphyCond(Float32.([3.4473, 3.6441, 3.8203, 3.9789, 4.1224, 4.2527,
            4.3714, 4.4800, 4.5796, 4.6711, 4.7554]),),
        Ni2011=RockphyCond(Float32.([-2.3579, -1.3975, -0.7500, -0.2847, 0.0652,
            0.3375, 0.5552, 0.7329, 0.8807, 1.0053, 1.1117]),),
        Sifre2014=RockphyCond(Float32.([-0.0128, 0.1570, 0.3119, 0.4537, 0.5837, 0.7031,
            0.8131, 0.9145, 1.0084, 1.0954, 1.1763]),),
        Gaillard2008=RockphyCond(Float32.([2.2276, 2.2577, 2.2865, 2.3140, 2.3403, 2.3655,
            2.3897, 2.4129, 2.4352, 2.4566, 2.4772]),))

    for i in eachindex(methods_list)
        m = keys(inps)[i]
        model = methods_list[i](inps[m]...)
        out_ = forward(model, [])
        @inferred forward(model, [])
        @test_opt forward(model, [])
        @test_call forward(model, [])
        for k in fieldnames(RockphyCond)
            @test all(isapprox.(getfield(out_, k), getfield(outs[m], k), rtol=1e-2))
        end
        modelD = methodsD_list[i](inps[m]...)
        @test sample_type(modelD) <: methods_list[i]
        params_ = default_params(typeof(model))
        @inferred forward(model, [], params_)
        @test_opt forward(model, [], params_)
        @test_call forward(model, [], params_)
    end
end
