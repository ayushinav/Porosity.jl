"""
    function inverse!(mₖ, robs, vars, alg_cache::occam_cache; 
            W, L, max_iters, χ2, response_fields, model_trans_utils,
            response_trans_utils, mᵣ, reg_term, verbose
        ):

updates `mₖ` using occam iteration to fit `robs` within a misfit of `χ2`, by default set to 1.0

## Arguments

  - `mₖ`: Initial model guess, will be updated during the inverse process
  - `robs`: response to invert for
  - `vars`: variables required for forward modeling, eg., `ω` for MT
  - `alg_cache`: deterimines the algorithm to be performed for inversion

## Keyword arguments

  - `W`: Weight matrix, defaults to identity matrix `I`
  - `L`: Regularization matrix, defaults to derivative matrix, given by `∂`(@ref)
  - `max_iters`: maximum number of iterations, defaults to 30
  - `χ2`: target misfit, defaults to 1.
  - `response_fields: choose data of response to perform inversion on, eg., ρₐ for MT, by default chooses all the data (ρₐ and ϕ)
  - `model_trans_utils:transform_utils`: conversion to and from computational domain, defaults to `log_sigmoid_tf`
  - `response_trans_utils`: for scaling the response parameters,
  - `mᵣ`: model in physical domain to be regularized against
  - `reg_term`: (For internals, only used by Occam) When model in physical domain does not exist, `reg_term` helps, eg. case of RTO-TKO
  - `verbose`: whether to print updates after each iteration, defaults to true

## Returns:

return message in the form of `return_code` and updates `mₖ` in-place.

## Usage:

```jldoctest
using LinearAlgebra
h = [1000.0, 1000.0] # m
ρ = log10.([100.0, 10.0, 1000.0]) # Ωm
m = MTModel(ρ, h)
T = 10 .^ (range(-2, 2; length=57))
ω = 2π ./ T
nω = length(T)
r_obs = forward(m, ω)
m_occam = MTModel(fill(2.0, 50), collect(range(0, 5e3; length=49)))
err_ρ = 0.1 .* r_obs.ρₐ
err_ϕ = asin(0.01) * 180 / π .* ones(length(ω))
err_resp = diagm(inv.([err_ρ..., err_ϕ...])) .^ 2
ret_code = inverse!(m_occam, r_obs, ω, Occam(; μgrid=[1e-2, 1e6]);
    W=err_resp, χ2=1.0, max_iters=50, verbose=true)

# output

1: Works golden section search: μ= 2559.998204059276, χ²= 22.867791875474534
2: Works golden section search: μ= 55.11264178720504, χ²= 16.89701371365216
3: Works golden section search: μ= 36.07944656829589, χ²= 12.58408431380094
4: Works golden section search: μ= 26.62727840194216, χ²= 9.43690878244996
5: Works golden section search: μ= 16.28674219250041, χ²= 7.1269845988502025
6: Works golden section search: μ= 8.309723956102522, χ²= 5.41891016151562
7: Works golden section search: μ= 4.116854777415263, χ²= 4.142014283721758
8: Works golden section search: μ= 2.115717196002266, χ²= 3.1786724887517708
9: Works golden section search: μ= 1.173804966446848, χ²= 2.447214296425664
10: Works golden section search: μ= 0.7087166896720183, χ²= 1.8893177576767075
11: Works golden section search: μ= 0.45988622588032063, χ²= 1.4623667912514455
12: Works golden section search: μ= 0.3163949884082966, χ²= 1.1347540808692684
13: Works golden section search: μ= 0.22874422627086552, χ²= 0.8828031760948966
MT.return_code{MTModel{Vector{Float64}, Vector{Float64}}}(true, (μ = 0.22874422627086552,), , 1.0, 0.8828031760948966)
```
"""
function inverse!(mₖ::model1,
        robs::response,
        vars::Vector{Float64},
        alg_cache::occam_cache;
        W=nothing,
        L=nothing,
        max_iters=30,
        χ2=1.0,
        response_fields::Vector{Symbol}=[k for k in fieldnames(typeof(robs))],
        model_trans_utils::trans_utils_T=sigmoid_tf,
        response_trans_utils::resp_utils_T=default_mt_tf_fns,
        mᵣ=nothing,
        reg_term=nothing,
        verbose::Union{Bool, Int}=true) where {model1 <: AbstractGeophyModel,
        response <: AbstractGeophyResponse, trans_utils_T, resp_utils_T}
    prec = eltype(mₖ.m)
    model_fields = [:m]

    n_model = length(mₖ.m)
    n_vars = length(vars)
    n_resp = length(response_fields) * n_vars

    (W === nothing) && (W = prec.(I(n_resp)))
    (L === nothing) && (L = prec.(∂(n_model)))

    respₖ = zero_abstract(robs)

    nresps = sum([length(getfield(robs, k)) for k in response_fields])
    nmods = length(mₖ.m)

    jc = zeros(eltype(mₖ.m), nresps, nmods)
    # jc = jacobian_cache(response_fields, robs, mₖ, model_fields).j

    lin_utils = linear_utils(view(mₖ.m, :), zeros(prec, n_resp), view(jc, :, :))

    for (i, k) in enumerate(response_fields)
        setfield!(respₖ, k, view(lin_utils.Fₖ, ((i - 1) * n_vars + 1):(i * n_vars)))
    end

    inv_utils = inverse_utils(
        L, W, reduce(vcat, [copy(getfield(robs, k)) for k in response_fields]))

    mₖ₊₁ = copy(mₖ)
    respₖ₊₁ = copy(respₖ)

    lin_prob = LinearProblem(inv_utils.D'inv_utils.D,
        lin_utils.Jₖ' * (inv_utils.dobs + lin_utils.Jₖ * lin_utils.mₖ - lin_utils.Fₖ))
    linsolve_prob = LinearSolve.init(lin_prob;
        assumptions=LinearSolve.OperatorAssumptions(
            true; condition=LinearSolve.OperatorCondition.WellConditioned))

    forward!(respₖ, mₖ, vars, response_trans_utils) # for the first iteration
    itr = 1
    chi2 = prec(1e6)

    if mᵣ !== nothing
        for k in model_fields # to computational domain
            getfield(mᵣ, k) .= model_trans_utils.itf.((getfield(mᵣ, k)))
        end
    end

    if reg_term === nothing
        reg_term = zero(mₖ.m)
    end

    μ_last = 0.0
    resp_cache = zero(robs)

    rvec = zero(lin_utils.Fₖ)

    model_type = typeof(mₖ).name.wrapper
    prep_j = prepare_jacobian(
        wrapper_DI!, rvec, AutoEnzyme(; mode=set_runtime_activity(Enzyme.Reverse)),
        mₖ.m, Constant(mₖ.h), Constant(vars), Cache(resp_cache),
        Constant(response_fields), Constant(response_trans_utils), Constant(model_type))

    DifferentiationInterface.jacobian!(wrapper_DI!, rvec, jc, prep_j,
        AutoEnzyme(; mode=set_runtime_activity(Enzyme.Reverse)), mₖ.m,
        Constant(mₖ.h), Constant(vars), Cache(resp_cache), Constant(response_fields),
        Constant(response_trans_utils), Constant(model_type))

    while itr <= max_iters
        do_verbose(itr, verbose) && (print("$itr: "))
        # @time MT.jacobian!(jc, mₖ, vars, model_fields, response_fields)
        # jc.j .= first(Enzyme.jacobian(set_runtime_activity(Reverse), f_temp, mₖ.m))

        DifferentiationInterface.jacobian!(
            wrapper_DI!, rvec, jc, AutoEnzyme(; mode=set_runtime_activity(Enzyme.Reverse)),
            mₖ.m, Constant(mₖ.h), Constant(vars),
            Cache(resp_cache), Constant(response_fields),
            Constant(response_trans_utils), Constant(model_type))

        for k in model_fields # to computational domain
            getfield(mₖ, k) .= model_trans_utils.itf.(getfield(mₖ, k))
        end

        μ_last = occam_step!(mₖ₊₁, # to store the next update, which will eventually be copied to mₖ
            respₖ₊₁, # to store the response for mₖ₊₁, for error calculation and anything
            vars, # to compute the forward model
            χ2, # threshold chi-squared error that needs to be met
            alg_cache.μgrid, # for gridsearch of μ for Occam
            lin_utils, # contains the mₖ, Jₖ, Fₖ associate with the current iteration
            inv_utils, # contains D= ∂(n), W and dobs
            model_trans_utils, # to  transform to and from the computational domain
            response_trans_utils, linsolve_prob; # for faster inverse operations
            model_fields=model_fields, response_fields=response_fields,
            mᵣ=mᵣ, verbose=verbose, reg_term=reg_term)

        mₖ.m .= mₖ₊₁.m

        forward!(respₖ, mₖ, vars, response_trans_utils)
        chi2 = χ²(reduce(vcat, [getfield(respₖ, k) for k in response_fields]),
            inv_utils.dobs; W=inv_utils.W)
        if chi2 < χ2
            break
        end
        itr += 1
    end

    if mᵣ !== nothing
        for k in model_fields # back to model domain
            getfield(mᵣ, k) .= model_trans_utils.tf.((getfield(mᵣ, k)))
        end
    end

    return return_code(chi2 <= χ2, (μ=μ_last,), mₖ, χ2, chi2)
end
