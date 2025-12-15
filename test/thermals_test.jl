@testitem "thermals : type inference" tags=[:thermals] begin
    using JET
    # solidus
    P = [3.0, 4.0, 5.0]
    ps_nt = (; P)
    @test_opt solidus_Hirschmann2000(ps_nt)
    @test_call solidus_Hirschmann2000(ps_nt)

    @test_opt solidus_Katz2003(ps_nt)
    @test_call solidus_Katz2003(ps_nt)

    # ΔT
    P = [3.0]
    T_solidus = [1400.0] .+ 273
    Ch2o_m = 0:200:1000
    Cco2_m = 0:200:1000
    ps_nt = (; P, T_solidus, Ch2o_m, Cco2_m)

    @test_opt ΔT_h2o_Katz2003(ps_nt)
    @test_call ΔT_h2o_Katz2003(ps_nt)

    @test_opt ΔT_co2_Dasgupta2007(ps_nt)
    @test_call ΔT_co2_Dasgupta2007(ps_nt)

    @test_opt ΔT_co2_Dasgupta2013(ps_nt)
    @test_call ΔT_co2_Dasgupta2013(ps_nt)

    @test_opt ΔT_h2o_Blatter2022(ps_nt)
    @test_call ΔT_h2o_Blatter2022(ps_nt)

    # volatile partition
    ϕ = 0.05
    Cco2 = 1000
    Cco2_sat = 38e4
    Ch2o = 1000
    ps_nt = (; ϕ, Cco2, Cco2_sat, Ch2o, D=0.0093)

    @test_opt get_Cco2_m(ps_nt)
    @test_call get_Cco2_m(ps_nt)

    @test_opt get_Ch2o_m(ps_nt)
    @test_call get_Ch2o_m(ps_nt)

    # melt fraction
    P = [3.0]
    T = (1440:10:1480) .+ 273.0
    T_solidus = [1450.0] .+ 273
    Ch2o = 100.0
    Cco2 = 1000.0

    ps_nt = (; P, T, T_solidus, Ch2o, Cco2)
    get_melt_fraction(ps_nt)
    @inferred get_melt_fraction(ps_nt)
    @test_call get_melt_fraction(ps_nt)
end
