"""
    rto_cache(m₀, μgrid, alg, max_iters, n_samples, χ2, response_fields, L, verbose)

returns `rto_cache` that specifies the algorithm to be used for stochastic inversion using RTO-TKO

### Arguments

  - `m₀`: `model`, usually a starting value to start inversion but mostly to allocate memory space
  - `μgrid`: prior space of μ
  - `alg`: type of algorithm for optimization step, choices include (`occam_cache`)[@ref], (`nl_cache`)[@ref], and (`opt_cache`)[@ref]
  - `max_iters`: maximum iterations for the optimization scheme
  - `n_samples`: number of samples to be obtained. Note: If there are any spurious samples (say `NaN` values), they will be automatically deleted, reducing the sample size
  - `χ2`: target misfit
  - `response_fields`: choose data of response to perform inversion on, eg., ρₐ for MT, by default chooses all the data (ρₐ and ϕ)
  - `L`: Regularization matrix, defaults to discretized derivative matrix given by ∂(@ref)
  - `verbose`: to print results or not
"""
struct rto_cache{model <: AbstractGeophyModel, T1 <: AbstractVector{<:Any},
    T2 <: Union{Float64, Float32}, A <: Union{occam_cache, nl_cache}, T3}
    m₀::model
    μgrid::T1
    alg::A
    max_iters::Int
    n_samples::Int
    χ2::T2
    response_fields::Vector{Symbol}
    L::T3
    verbose::Bool
end

function rto_cache(m₀, μgrid, alg, max_iters, n_samples, χ2, response_fields, verbose)
    return rto_cache(
        m₀, μgrid, alg, max_iters, n_samples, χ2, response_fields, ∂(length(m₀.m)), verbose)
end

"""
    stochastic_inverse(r_obs::resp1,
            err_resp::resp2,
            vars,
            alg_cache::rto_cache;
            model_trans_utils::NamedTuple=(m=sigmoid_tf, h=lin_tf),
            response_trans_utils::NamedTuple=(ρₐ=lin_tf, ϕ=lin_tf)) where {
            resp1 <: AbstractGeophyResponse, resp2 <: AbstractGeophyResponse}

### Returns

    `AbstractMCMC.jl Chain` containing the vectors used for performing samping. Note that all the variables will be named `var_inf`, for variable inferred. 
    The vector in each sample will be all the fields of the model being inferred on concatenated together.

### Arguments

  - `r_obs`: `response` that needs to inverted for
  - `err_resp`: `response` variable containing the errors associated with observed response
  - `vars`: variables that need to be passed into the `forward` function along with `model` to generate a `response`
  - `alg_cache`: to tell the compiler what type of stochastic inversion method is to be used

### Keyword Arguments

  - `model_trans_utils`: A named tuple containing `transform_utils` for the fields of model that need to be scaled/modified. If not provided for any `model` field, the field won't be modified
  - `response_trans_utils`: for scaling the response parameters
"""
function stochastic_inverse(r_obs::resp1,
        err_resp::resp2,
        vars,
        alg_cache::rto_cache;
        model_trans_utils::NamedTuple=(m=sigmoid_tf, h=lin_tf),
        response_trans_utils::NamedTuple=(ρₐ=lin_tf, ϕ=lin_tf),
        progress_bar=true) where {
        resp1 <: AbstractGeophyResponse, resp2 <: AbstractGeophyResponse}
    W = Diagonal(vcat([inv.(getfield(err_resp, k)) .^ 2
                       for k in fieldnames(typeof(err_resp))]...))

    reg_term = copy(alg_cache.m₀.m)
    pert_resp = copy(r_obs)

    n = length(alg_cache.m₀.m)
    L = ∂(n)

    μ = 1.0
    resp_ = copy(r_obs)
    ϕ = (1 + sqrt(5)) / 2

    m_chains = zeros(n, alg_cache.n_samples)
    μ_chains = zeros(1, alg_cache.n_samples)

    i = 1
    # (!(verbose==false)) && (prog = Progress(alg_cache.n_samples; enabled=true))
    (progress_bar) && (prog = Progress(alg_cache.n_samples; enabled=true))

    while i <= (alg_cache.n_samples)

        ## Step 1

        # perturbed response
        for k in fieldnames(typeof(r_obs))
            getfield(pert_resp, k) .= getfield(r_obs, k) .+
                                      randn(size(getfield(err_resp, k))) .*
                                      getfield(err_resp, k)
        end

        # perturbed model : the one we regularize against
        # we first draw a perturbation in the computational domain

        #reg term
        mul!(reg_term, L'L, randn(eltype(reg_term), size(reg_term))) # L'L * ξ
        rmul!(reg_term, μ) # μ * L'L * ξ 

        fill!(alg_cache.m₀.m, 2.0)
        # @show norm(alg_cache.m₀.m)

        if typeof(alg_cache.alg) <: occam_cache
            ret_code = inverse!(alg_cache.m₀, pert_resp, vars, Occam(; μgrid=[μ, μ]);
                W=W, χ2=alg_cache.χ2, max_iters=alg_cache.max_iters,
                response_fields=alg_cache.response_fields, verbose=alg_cache.verbose,
                reg_term=reg_term, model_trans_utils=model_trans_utils[:m],
                response_trans_utils=response_trans_utils)

        elseif typeof(alg_cache.alg) <: nl_cache
            ret_code = inverse!(
                alg_cache.m₀, pert_resp, vars, NonlinearAlg(; alg=alg_cache.alg.alg, μ=μ);
                W=W, χ2=alg_cache.χ2, L=alg_cache.L, max_iters=alg_cache.max_iters,
                response_fields=alg_cache.response_fields,
                verbose=alg_cache.verbose, model_trans_utils=model_trans_utils[:m])
        elseif typeof(alg_cache.alg) <: opt_cache
            ret_code = inverse!(
                alg_cache.m₀, pert_resp, vars, OptAlg(; alg=alg_cache.alg.alg, μ=μ);
                W=W, χ2=alg_cache.χ2, L=alg_cache.L, max_iters=alg_cache.max_iters,
                response_fields=alg_cache.response_fields,
                verbose=alg_cache.verbose, model_trans_utils=model_trans_utils[:m])
        end

        m_chains[:, i] .= alg_cache.m₀.m

        ## Step 2

        # perturbed response
        for k in fieldnames(typeof(r_obs))
            getfield(pert_resp, k) .= getfield(r_obs, k) .+
                                      randn(size(getfield(err_resp, k))) .*
                                      getfield(err_resp, k)
        end

        # broadcast!((x) -> (model_trans_utils[:m].itf(x)), alg_cache.m₀.m, alg_cache.m₀.m) # to computational domain
        rmul!(alg_cache.m₀.m, sqrt(μ)) # current μ

        function g(x) #, alg_cache, model_trans_utils, resp_, vars)
            rmul!(alg_cache.m₀.m, sqrt(inv(x)))

            # broadcast!(model_trans_utils[:m].itf, alg_cache.m₀.m, alg_cache.m₀.m) # to model domain
            # broadcast!((x) -> (10.0^model_trans_utils[:m].tf(x)), alg_cache.m₀.m, alg_cache.m₀.m) #  to model domain
            # alg_cache.m₀.m .= model_trans_utils[:m].tf.(alg_cache.m₀.m) # to model domain

            forward!(resp_, alg_cache.m₀, vars, response_trans_utils)
            # @show "chi2_err"
            chi2_err = χ²(
                reduce(vcat, [getfield(resp_, k) for k in alg_cache.response_fields]),
                reduce(vcat, [getfield(pert_resp, k) for k in alg_cache.response_fields]);
                W=W)

            # broadcast!(trans_utils[:m].tf, alg_cache.m₀.m, alg_cache.m₀.m) # to computational domain
            # alg_cache.m₀.m .= model_trans_utils[:m].itf.(alg_cache.m₀.m) # to computational domain
            # broadcast!(
            #     (x) -> (trans_utils[:m].itf(log10.(x))), alg_cache.m₀.m, alg_cache.m₀.m) # to computational domain
            rmul!(alg_cache.m₀.m, sqrt(x))

            return chi2_err
        end

        x₁ = alg_cache.μgrid[1]
        x₃ = alg_cache.μgrid[end]
        x₂ = 10.0^((log10(x₃) + ϕ * log10(x₁)) / (1 + ϕ))
        x₄ = 10.0^((log10(x₁) + ϕ * log10(x₃)) / (1 + ϕ))

        fx₂ = g(x₂)
        fx₄ = g(x₄)

        tol = 1e-5
        count = 0
        while (x₃ - x₁) >= tol
            count += 1
            if count > 100
                print("100 golden section iterations done. \t")
                break
            end
            if fx₄ > fx₂
                x₃ = x₄
                x₄ = x₂
                fx₄ = fx₂
                x₂ = 10.0^((log10(x₃) + ϕ * log10(x₁)) / (1 + ϕ))
                fx₂ = g(x₂)

            else
                x₁ = x₂
                x₂ = x₄
                fx₂ = fx₄
                x₄ = 10.0^((log10(x₁) + ϕ * log10(x₃)) / (1 + ϕ))
                fx₄ = g(x₄)
            end
        end
        μ = sqrt(x₁ * x₃)

        μ_chains[1, i] = μ

        i += 1

        (progress_bar) && (next!(prog; showvalues=[(Symbol("#samples"), i)]))
        # @show i
    end

    idcs = broadcast(!isnan, view(m_chains, 1, :))
    @show sum(idcs)

    return Turing.Chains(
        vcat(m_chains[:, idcs], μ_chains[:, idcs])', [Symbol("m[$i]") for i in 1:(n + 1)])
end

# mutable struct RTO_MTModel <: AbstractGeophyModel
#     ξ
#     h
#     μ
# end
