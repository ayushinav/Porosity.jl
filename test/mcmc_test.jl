@testitem "fixed discretization : NUTS" tags=[:mcmc] begin
    using Distributions, Turing, LinearAlgebra

    m_test = MTModel(log10.([100.0, 10.0, 1000.0]), [1e3, 1e3])
    f = 10 .^ range(-4; stop=1, length=25)
    ω = vec(2π .* f)

    r_obs = forward(m_test, ω)

    err_phi = asin(0.01) * 180 / π .* ones(length(ω))
    err_appres = 0.02 * r_obs.ρₐ
    err_resp = MTResponse(err_appres, err_phi)

    r_obs.ρₐ .= r_obs.ρₐ .+ err_appres
    r_obs.ϕ .= r_obs.ϕ .+ err_phi

    respD = MTResponseDistribution(normal_dist, normal_dist)

    z = 10 .^ collect(range(1; stop=4, length=50))
    h = diff(z)

    modelD = MTModelDistribution(
        product_distribution([Uniform(-1.0, 5.0) for i in eachindex(z)]),
        vec(h)
    )

    n_samples = 20
    mcache = mcmc_cache(modelD, respD, n_samples, NUTS())

    mcmc_chain = stochastic_inverse(r_obs, err_resp, ω, mcache)

    model_list = get_model_list(mcmc_chain, modelD)
    W = diagm(inv.([err_resp.ρₐ..., err_resp.ϕ...])) .^ 2

    err = zeros(n_samples)
    for idx in 1:n_samples
        m_model = model_list[idx]
        resp_model = forward(m_model, ω)

        err[idx] = χ²(
            reduce(vcat, [getfield(resp_model, k) for k in [:ρₐ, :ϕ]]),
            reduce(vcat, [getfield(r_obs, k) for k in [:ρₐ, :ϕ]]);
            W=W
        )
    end
    @show err

    @test sum(err[11:20]) <= sum(err[1:10])
end

@testitem "variable discretization : NUTS" tags=[:mcmc] begin
    using Distributions, Turing, LinearAlgebra

    m_test = MTModel(log10.([100.0, 10.0, 1000.0]), [1e3, 1e3])
    f = 10 .^ range(-4; stop=1, length=25)
    ω = vec(2π .* f)

    r_obs = forward(m_test, ω)

    err_phi = asin(0.01) * 180 / π .* ones(length(ω))
    err_appres = 0.02 * r_obs.ρₐ
    err_resp = MTResponse(err_appres, err_phi)

    r_obs.ρₐ .= r_obs.ρₐ .+ err_appres
    r_obs.ϕ .= r_obs.ϕ .+ err_phi

    respD = MTResponseDistribution(normal_dist, normal_dist)

    z = 10 .^ collect(range(1; stop=4, length=50))
    h = diff(z)
    h_bounds = [[ih / 3, ih * 3] for ih in h]

    modelD = MTModelDistribution(
        product_distribution([Uniform(-1.0, 5.0) for i in eachindex(z)]),
        product_distribution([Uniform(ih_bounds...) for ih_bounds in h_bounds])
    )

    n_samples = 20
    mcache = mcmc_cache(modelD, respD, n_samples, NUTS())

    mcmc_chain = stochastic_inverse(r_obs, err_resp, ω, mcache)

    model_list = get_model_list(mcmc_chain, modelD)
    W = diagm(inv.([err_resp.ρₐ..., err_resp.ϕ...])) .^ 2

    err = zeros(n_samples)
    for idx in 1:n_samples
        m_model = model_list[idx]
        resp_model = forward(m_model, ω)

        err[idx] = χ²(
            reduce(vcat, [getfield(resp_model, k) for k in [:ρₐ, :ϕ]]),
            reduce(vcat, [getfield(r_obs, k) for k in [:ρₐ, :ϕ]]);
            W=W
        )
    end
    @show err

    @test sum(err[11:20]) <= sum(err[1:10])
end

@testitem "fixed discretization : SliceSampler" tags=[:mcmc] begin
    using Distributions, Pigeons, LinearAlgebra

    m_test = MTModel(log10.([100.0, 10.0, 1000.0]), [1e3, 1e3])
    f = 10 .^ range(-4; stop=1, length=25)
    ω = vec(2π .* f)

    r_obs = forward(m_test, ω)

    err_phi = asin(0.01) * 180 / π .* ones(length(ω))
    err_appres = 0.02 * r_obs.ρₐ
    err_resp = MTResponse(err_appres, err_phi)

    r_obs.ρₐ .= r_obs.ρₐ .+ err_appres
    r_obs.ϕ .= r_obs.ϕ .+ err_phi

    respD = MTResponseDistribution(normal_dist, normal_dist)

    z = 10 .^ collect(range(1; stop=4, length=50))
    h = diff(z)

    modelD = MTModelDistribution(
        product_distribution([Uniform(-1.0, 5.0) for i in eachindex(z)]),
        vec(h)
    )

    n_samples = 16
    mcache = mcmc_cache(modelD, respD, n_samples, SliceSampler())

    mcmc_chain = stochastic_inverse(r_obs, err_resp, ω, mcache)

    model_list = get_model_list(mcmc_chain, modelD)
    W = diagm(inv.([err_resp.ρₐ..., err_resp.ϕ...])) .^ 2

    err = zeros(n_samples)
    for idx in 1:n_samples
        m_model = model_list[idx]
        resp_model = forward(m_model, ω)

        err[idx] = χ²(
            reduce(vcat, [getfield(resp_model, k) for k in [:ρₐ, :ϕ]]),
            reduce(vcat, [getfield(r_obs, k) for k in [:ρₐ, :ϕ]]);
            W=W
        )
    end
    @show err

    @test sum(err[9:16]) <= sum(err[1:8])
end

@testitem "variable discretization  :SliceSampler" tags=[:mcmc] begin
    using Distributions, Pigeons, LinearAlgebra

    m_test = MTModel(log10.([100.0, 10.0, 1000.0]), [1e3, 1e3])
    f = 10 .^ range(-4; stop=1, length=25)
    ω = vec(2π .* f)

    r_obs = forward(m_test, ω)

    err_phi = asin(0.01) * 180 / π .* ones(length(ω))
    err_appres = 0.02 * r_obs.ρₐ
    err_resp = MTResponse(err_appres, err_phi)

    r_obs.ρₐ .= r_obs.ρₐ .+ err_appres
    r_obs.ϕ .= r_obs.ϕ .+ err_phi

    respD = MTResponseDistribution(normal_dist, normal_dist)

    z = 10 .^ collect(range(1; stop=4, length=50))
    h = diff(z)
    h_bounds = [[ih / 3, ih * 3] for ih in h]

    modelD = MTModelDistribution(
        product_distribution([Uniform(-1.0, 5.0) for i in eachindex(z)]),
        product_distribution([Uniform(ih_bounds...) for ih_bounds in h_bounds])
    )

    n_samples = 16
    mcache = mcmc_cache(modelD, respD, n_samples, SliceSampler())

    mcmc_chain = stochastic_inverse(r_obs, err_resp, ω, mcache)

    model_list = get_model_list(mcmc_chain, modelD)
    W = diagm(inv.([err_resp.ρₐ..., err_resp.ϕ...])) .^ 2

    err = zeros(n_samples)
    for idx in 1:n_samples
        m_model = model_list[idx]
        resp_model = forward(m_model, ω)

        err[idx] = χ²(
            reduce(vcat, [getfield(resp_model, k) for k in [:ρₐ, :ϕ]]),
            reduce(vcat, [getfield(r_obs, k) for k in [:ρₐ, :ϕ]]);
            W=W
        )
    end
    @show err

    @test sum(err[9:16]) <= sum(err[1:8])
end
