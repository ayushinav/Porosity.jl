"""
`occam_cache`: specifies the inverse algorithm while having a cache.
"""
mutable struct occam_cache{T}
    μgrid::Vector{T}
end
"""
`linsolve!`: Performs `inv(B)*y` using `LinearSolve.jl`
"""
function linsolve!(x, prob_init, B, y)
    prob_init.A = B
    prob_init.b = y
    x .= solve!(prob_init)
    nothing
end
"""
`Occam(;μgrid= [0.01, 1e6])`
"""
function Occam(; μgrid=[0.01, 1e6])
    return occam_cache{eltype(μgrid)}(μgrid)
end

"""
    function occam_step!(mₖ₊₁::model,
        respₖ₊₁::response,
        vars::Union{AbstractVector{Float32}, AbstractVector{Float64}},
        χ2::Union{Float64, Float32},
        μgrid::Vector{Float64},
        lin_utils::linear_utils,
        inv_utils::inverse_utils,
        model_trans_utils::transform_utils,
        response_trans_utils::NamedTuple,
        linsolve_prob::LinearSolve.LinearCache;
        model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(mₖ₊₁))],
        response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(respₖ₊₁))],
        verbose= false
        ):

performs a single step of occam inversion, using golden line search.

## Arguments

  - `mₖ`: Initial model guess, will be updated during the inverse process
  - `respₖ₊₁`: response to invert for
  - `vars`: variables required for forward modeling, eg., `ω` for MT
  - `alg_cache`: deterimines the algorithm to be performed for inversion
  - `W`: Weight matrix, defaults to identity matrix `I`
  - `L`: Regularization matrix, defaults to discretized derivative matrix given by ∂(@ref)
  - `max_iters`: maximum number of iterations, defaults to 30
  - `χ2`: target misfit, defaults to 1.0
  - `response_fields: choose data of response to perform inversion on, eg., ρₐ for MT, by default chooses all the data (ρₐ and ϕ)
  - `model_trans_utils`: conversion to and from computational domain,
  - `response_trans_utils`: NamedTuple containing `transform_utils` for scaling different response parameters,
  - `mᵣ`: model in physical domain to be regularized against
  - `reg_term`: (For internals) When model in physical domain does not exist, `reg_term` helps, eg. case of RTO-TKO
  - `verbose`: whether to print updates after each iteration, defaults to true
"""
function occam_step!(mₖ₊₁::model1, # to store the next update, which will eventually be copied to mₖ
        respₖ₊₁::response, # to store the response for mₖ₊₁, for error calculation and anything
        vars::Union{AbstractVector{Float32}, AbstractVector{Float64}}, # to compute the forward model
        χ2::Union{Float64, Float32}, # threshold chi-squared error that needs to be met
        μgrid::Vector{Float64}, # contains end points of the bounds for the lagrange multiplier
        lin_utils::linear_utils, # contains the mₖ, Jₖ, Fₖ associate with the current iteration
        inv_utils::inverse_utils, # contains D= ∂(n), W and dobs
        model_trans_utils::transform_utils, # to  transform to and from the computational domain
        response_trans_utils::NamedTuple, #for scaling the response parameters,
        linsolve_prob::LinearSolve.LinearCache; # for faster inverse operations
        model_fields::Vector{Symbol}=[k for k in fieldnames(typeof(mₖ₊₁))],
        response_fields::Vector{Symbol}=[k for k in fieldnames(typeof(respₖ₊₁))],
        mᵣ::model2,
        verbose::Bool=true,
        reg_term::AbstractVector) where {
        model1 <: AbstractGeophyModel, model2 <: Union{AbstractGeophyModel, Nothing},
        response <: AbstractGeophyResponse}
    ϕ = (1 + sqrt(5)) / 2
    chi2min = (typeof(χ2))(1e6)
    μ = zero(eltype(μgrid))
    count = 0 # so that iterations do not run forever (will rarely happen, if it will)

    function f(x, mᵣ)
        find_x(x, mₖ₊₁, respₖ₊₁, vars, inv_utils, lin_utils, model_fields, response_fields,
            model_trans_utils, response_trans_utils, linsolve_prob, reg_term, mᵣ)
    end

    x₁ = μgrid[1]
    x₃ = μgrid[end]
    x₂ = 10.0^((log10(x₃) + ϕ * log10(x₁)) / (1 + ϕ))
    x₄ = 10.0^((log10(x₁) + ϕ * log10(x₃)) / (1 + ϕ))

    # fx₁ = f(x₁)
    # fx₃ = f(x₃)
    fx₂ = f(x₂, mᵣ)
    fx₄ = f(x₄, mᵣ)

    tol = 1e-5
    count = 0
    while (x₃ - x₁) >= tol
        count += 1
        if count > 100
            verbose && (print("100 golden section iterations done. \t"))
            break
        end
        if fx₄ > fx₂
            x₃ = x₄
            x₄ = x₂
            fx₄ = fx₂
            x₂ = 10.0^((log10(x₃) + ϕ * log10(x₁)) / (1 + ϕ))
            fx₂ = f(x₂, mᵣ)

        else
            x₁ = x₂
            x₂ = x₄
            fx₂ = fx₄
            x₄ = 10.0^((log10(x₁) + ϕ * log10(x₃)) / (1 + ϕ))
            fx₄ = f(x₄, mᵣ)
        end
    end
    μ = sqrt(x₁ * x₃)

    # At the moment mₖ₊₁ contains the update for the last μ, we rewrite it with the best μ found.

    if mᵣ === nothing
        linsolve!(mₖ₊₁.m,
            linsolve_prob,
            μ .* inv_utils.D' * inv_utils.D .+ lin_utils.Jₖ' * inv_utils.W * lin_utils.Jₖ,
            lin_utils.Jₖ' *
            inv_utils.W *
            (inv_utils.dobs + lin_utils.Jₖ * lin_utils.mₖ - lin_utils.Fₖ) + reg_term)
    else
        linsolve!(mₖ₊₁.m,
            linsolve_prob,
            μ .* inv_utils.D' * inv_utils.D .+ lin_utils.Jₖ' * inv_utils.W * lin_utils.Jₖ,
            lin_utils.Jₖ' *
            inv_utils.W *
            (inv_utils.dobs + lin_utils.Jₖ * lin_utils.mₖ - lin_utils.Fₖ) +
            μ .* inv_utils.D' * inv_utils.D * mᵣ.m +
            reg_term)
    end

    for k in model_fields # to model domain
        getfield(mₖ₊₁, k) .= model_trans_utils.tf.(getfield(mₖ₊₁, k))
    end

    forward!(respₖ₊₁, mₖ₊₁, vars, response_trans_utils)

    do_verbose(verbose) && (print("Works golden section search: μ= $μ, χ²= ",
        χ²(reduce(vcat, [copy(getfield(respₖ₊₁, k)) for k in response_fields]),
            inv_utils.dobs; W=inv_utils.W),
        "\n"))
    return μ
end

function find_x(x::T1, mₖ₊₁::model, respₖ₊₁::response, vars, inv_utils::MT.inverse_utils,
        lin_utils::MT.linear_utils, model_fields::Vector{Symbol},
        response_fields::Vector{Symbol}, model_trans_utils::T3, response_trans_utils::T,
        linsolve_prob, reg_term, mᵣ::Nothing) where {T1, T, T3, model, response}
    linsolve!(mₖ₊₁.m,
        linsolve_prob,
        x .* inv_utils.D' * inv_utils.D .+ lin_utils.Jₖ' * inv_utils.W * lin_utils.Jₖ,
        lin_utils.Jₖ' *
        inv_utils.W *
        (inv_utils.dobs + lin_utils.Jₖ * lin_utils.mₖ - lin_utils.Fₖ) .+ reg_term)

    broadcast!(model_trans_utils.tf, mₖ₊₁.m, mₖ₊₁.m)
    forward!(respₖ₊₁, mₖ₊₁, vars, response_trans_utils)

    return χ²(reduce(vcat, [getfield(respₖ₊₁, k) for k in response_fields]),
        inv_utils.dobs; W=inv_utils.W)
end

function find_x(x::T1, mₖ₊₁::model, respₖ₊₁::response, vars, inv_utils::MT.inverse_utils,
        lin_utils::MT.linear_utils, model_fields::Vector{Symbol},
        response_fields::Vector{Symbol}, model_trans_utils::T3, response_trans_utils::T,
        linsolve_prob, reg_term, mᵣ) where {T1, T, T3, model, response}
    linsolve!(mₖ₊₁.m,
        linsolve_prob,
        x .* inv_utils.D' * inv_utils.D .+ lin_utils.Jₖ' * inv_utils.W * lin_utils.Jₖ,
        lin_utils.Jₖ' *
        inv_utils.W *
        (inv_utils.dobs + lin_utils.Jₖ * lin_utils.mₖ - lin_utils.Fₖ) +
        μ .* inv_utils.D' * inv_utils.D * mᵣ.m .+ reg_term)

    broadcast!(model_trans_utils.tf, mₖ₊₁.m, mₖ₊₁.m)
    forward!(respₖ₊₁, mₖ₊₁, vars, response_trans_utils)

    return χ²(reduce(vcat, [getfield(respₖ₊₁, k) for k in response_fields]),
        inv_utils.dobs; W=inv_utils.W)
end
