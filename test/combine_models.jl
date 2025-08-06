@testitem "combine models" tags = [:combine_models] begin
    using JET
    m = multi_rp_modelType(SEO3, anharmonic, Nothing, Nothing)
    ps_nt = (;
        T = [800.0f0, 1000.0f0] .+ 273,
        P = 3.0f0,
        ρ = 3300.0f0,
        Ch2o_m = 1000.0f0,
        ϕ = 0.1f0,
    )
    model = m(ps_nt)
    resp = to_resp_nt(forward(model, []))
    @inferred forward(model, [])

    resp_SEO3 = forward(model.cond, [])
    resp_anharmonic = forward(model.elastic, [])

    resp_check = merge(to_nt.([resp_SEO3, resp_anharmonic])...)

    @test resp == resp_check
end
