"""
`nl_cache`: specifies the inverse algorithm while having a cache.
"""
mutable struct nl_cache{T1, T2}
    alg::T1
    μ::T2
end
"""
    NonlinearAlg(; alg = LevenbergMarquardt, μ = 1.0)

returns `nl_cache` that specifies which non linear solver to use for the inverse problem

## Keyword Arguments

  - `alg`: `NonlinearSolve`[@ref] algorithm to be used, defaults to LevenbergMarquardt
  - `μ` : regularization weight
"""
function NonlinearAlg(; alg=LevenbergMarquardt, μ=1.0)
    return nl_cache(alg, μ)
end

# ===== NonlinearSolve.jl =========

function inverse!(mₖ::model1,
        robs::response,
        vars::Vector{Float64},
        alg_cache::nl_cache;
        W=nothing,
        L=nothing,
        max_iters=30,
        χ2=1.0,
        response_fields::Vector{Symbol}=[k for k in fieldnames(typeof(robs))],
        model_trans_utils::transform_utils=sigmoid_tf,
        response_trans_utils::NamedTuple=default_mt_tf_fns,
        mᵣ=nothing,
        reg_term=nothing,
        verbose::Union{Bool, Int}=true) where {
        model1 <: AbstractGeophyModel, response <: AbstractGeophyResponse}
    prec = eltype(mₖ.m)

    n_vars = length(vars)
    n_resp = length(response_fields) * n_vars

    (W === nothing) && (W = prec.(I(n_resp)))
    (L === nothing) && (L = prec.(∂(length(mₖ.m))))

    model_type = typeof(mₖ).name.wrapper
    (mᵣ === nothing) && (mᵣ = model_type(zero(mₖ.m), mₖ.h))

    p = (model_type=model_type, h=mₖ.h, model_trans_utils=model_trans_utils,
        response_trans_utils=response_trans_utils, vars=vars,
        response_fields=response_fields, W=W, μ=alg_cache.μ, r_obs=robs, L=L, mᵣ=mᵣ)

    prob = SciMLBase.NonlinearLeastSquaresProblem(
        construct_cost_function_for_nl_inv, model_trans_utils.itf.(mₖ.m), p)
    nlcache = init(prob, alg_cache.alg())
    iters = 1
    chi2 = 1e6

    resp_ = zero(robs)

    while iters <= max_iters
        model = model_type(model_trans_utils.tf.(nlcache.u), mₖ.h)
        forward!(resp_, model, vars, response_trans_utils)
        chi2 = χ²(reduce(vcat, [getfield(resp_, k) for k in response_fields]),
            reduce(vcat, [getfield(robs, k) for k in response_fields]); W=W)
        do_verbose(iters, verbose) && println("iteration = $iters : data misfit => $chi2")

        if chi2 <= χ2 # check misfit condition
            break
        end

        step!(nlcache)
        mₖ.m .= model_trans_utils.tf.(nlcache.u)
        iters += 1
    end

    return return_code(chi2 <= χ2, (μ=alg_cache.μ,), mₖ, χ2, chi2)
end

# focussing on just geophysical models for now
# Not performant at the moment
function construct_cost_function_for_nl_inv(m, p)
    @unpack model_type, h, model_trans_utils, response_trans_utils, vars, response_fields, W, μ, r_obs, L, mᵣ = p
    # model = model_type(model_trans_utils.tf.(m), h)
    model = model_type(broadcast(model_trans_utils.tf, m), h)
    resp_ = forward(model, vars, response_trans_utils)

    L1 = χ²(reduce(vcat, [getfield(resp_, k) for k in response_fields]),
        reduce(vcat, [getfield(r_obs, k) for k in response_fields]); W=W) * sqrt(size(W, 1))
    L2 = μ * norm(L * (model.m .- mᵣ.m))
    # @show L1, L2

    return [L1^2 + L2]
end
