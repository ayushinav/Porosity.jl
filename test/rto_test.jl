@testitem "RTO" tags=[:rto] begin
    using MT
    using Distributions, Turing, LinearAlgebra

    m_test = MTModel(log10.([100.0, 10.0, 1000.0]), [1e3, 1e3])
    f = 10 .^ range(-2; stop=2, length=57)
    ω = vec(2π .* f)

    r_obs = forward(m_test, ω)

    err_phi = asin(0.01) * 180 / π .* ones(length(ω))
    err_appres = 0.1 * r_obs.ρₐ
    err_resp = MTResponse(err_appres, err_phi)

    r_obs.ρₐ .= r_obs.ρₐ .+ rand(length(ω)) .* err_appres
    r_obs.ϕ .= r_obs.ϕ .+ rand(length(ω)) .* err_phi

    respD = MTResponseDistribution(normal_dist, normal_dist)

    z = collect(range(0, 5e3; length=50))
    h = diff(z)

    modelD = MTModelDistribution(
        product_distribution([Uniform(-1.0, 5.0) for i in eachindex(z)]),
        vec(h)
    )

    n_samples = 10
    m_rto = MTModel(2 .* ones(length(z)), vec(h))
    r_cache = rto_cache(m_rto, [1e-6, 1e6], Occam(), 50, n_samples, 1.0, [:ρₐ, :ϕ], false)

    rto_chain = stochastic_inverse(
        r_obs,
        err_resp,
        ω,
        r_cache;
        model_trans_utils=(; m=MT.lin_tf),
        progress_bar=true
    )

    mt_chain = Turing.Chains(
        (rto_chain.value.data[:, 1:length(z), :]),
        [Symbol("ρ[$i]") for i in 1:length(z)]
    )

    model_list = get_model_list(mt_chain, modelD)

    m_model = model_list[end]
    resp_model = forward(m_model, ω)

    W = diagm(inv.([err_appres..., err_phi...])) .^ 2

    err = sqrt(
        norm(
        inv(2 * length(ω)) *
        ([resp_model.ρₐ..., resp_model.ϕ...] .- [r_obs.ρₐ..., r_obs.ϕ...]) ./
        [err_resp.ρₐ..., err_resp.ϕ...],
        2
    ),
    )

    @test err <= 2.0
end
