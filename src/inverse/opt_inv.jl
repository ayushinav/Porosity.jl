"""
`nl_cache`: specifies the inverse algorithm while having a cache.
"""
mutable struct opt_cache{T1, T2}
    alg::T1
    μ::T2
end
"""
    OptAlg(; alg = LBFGS, μ = 1.0, kwargs...)

returns `opt_cache` that specifies which `Optimization.jl` solver to use for the inverse problem

## Keyword Arguments

  - `alg`: `Optimization`[@ref] algorithm to be used, defaults to LBFGS
  - `μ` : regularization weight
"""
function OptAlg(; alg=LBFGS, μ=1.0)
    return opt_cache(alg, μ)
end
# ======================== using Optimization.jl ===============================

function inverse!(mₖ::model1,
        robs::response,
        vars::Vector{Float64},
        alg_cache::opt_cache;
        W=nothing,
        L=nothing,
        max_iters=30,
        χ2=1.0,
        response_fields::Vector{Symbol}=[k for k in fieldnames(typeof(robs))],
        model_trans_utils::transform_utils=sigmoid_tf,
        response_trans_utils::NamedTuple=(; ρₐ=lin_tf, ϕ=lin_tf),
        mᵣ=nothing,
        reg_term=nothing,
        verbose::Union{Bool, Int}=true) where {
        model1 <: AbstractGeophyModel, response <: AbstractGeophyResponse}
    prec = eltype(mₖ.m)
    model_fields = [:m]

    n_vars = length(vars)
    n_resp = length(response_fields) * n_vars

    (W === nothing) && (W = prec.(I(n_resp)))
    (L === nothing) && (L = prec.(∂(length(mₖ.m))))

    model_type = typeof(mₖ).name.wrapper
    (mᵣ === nothing) && (mᵣ = model_type(zero(mₖ.m), mₖ.h))

    p = (model_type=model_type, h=mₖ.h, model_trans_utils=model_trans_utils,
        response_trans_utils=response_trans_utils, vars=vars,
        response_fields=response_fields, W=W, μ=alg_cache.μ, r_obs=robs, L=L, mᵣ=mᵣ)

    optfn = OptimizationFunction(
        construct_cost_function_for_opt, Optimization.AutoForwardDiff())
    prob = OptimizationProblem(optfn, model_trans_utils.itf.(mₖ.m), p)

    cb(state, l) = cb_(state, l, verbose, L, alg_cache.μ, model_trans_utils, χ2)
    sol = solve(prob, alg_cache.alg(); callback=cb, maxiters=max_iters)

    mₖ.m .= model_trans_utils.tf.(sol.u)

    resp_ = forward(mₖ, vars, response_trans_utils)
    chi2 = χ²(reduce(vcat, [getfield(resp_, k) for k in response_fields]),
        reduce(vcat, [getfield(robs, k) for k in response_fields]); W=W)

    return return_code(chi2 <= χ2, (μ=alg_cache.μ,), mₖ, χ2, chi2)
end

function cb_(state, l, verbose, L, μ, model_trans_utils, χ2)
    chi2 = sqrt(l - μ * norm(L * model_trans_utils.tf.(state.u)))
    do_verbose(verbose) && println("iteration = $(state.iter) : data misfit => $chi2")

    return (chi2 < χ2)
end

function construct_cost_function_for_opt(m, p)
    @unpack model_type, h, model_trans_utils, response_trans_utils, vars, response_fields, W, μ, r_obs, L, mᵣ = p
    # model = model_type(model_trans_utils.tf.(m), h)
    model = model_type(broadcast(model_trans_utils.tf, m), h)
    resp_ = forward(model, vars, response_trans_utils)

    L1 = χ²(reduce(vcat, [getfield(resp_, k) for k in response_fields]),
        reduce(vcat, [getfield(r_obs, k) for k in response_fields]); W=W) * sqrt(size(W, 1))
    L2 = μ * norm(L * (model.m .- mᵣ.m))
    # @show L1, L2

    return L1^2 + L2
end
