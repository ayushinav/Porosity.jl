"""
    stochastic_inverse(
        r_obs::response,
        err_resp::response,
        vars,
        alg_cache::mcmc_cache;
        model_trans_utils::NamedTuple = (m = lin_tf, h = lin_tf)
        )

function to perform sampling

## Returns

    `Chain` containing the samples. Note that all the variables will be named `m0`. If 

## Arguments

  - `r_obs`: `response` that needs to inverted for
  - `err_resp`: `response` variable containing the errors associated with observed response
  - `vars`: variables that need to be passed into the `forward` function along with `model` to generate a `response`
  - `alg_cache`: to tell the compiler what type of stochastic inversion method is to be used

### Optional Keyword Arguments

  - `n_chains` : Number of chains to use
  - `model_trans_utils`: A named tuple containing `transform_utils` for the fields of model that need to be scaled/modified,
    defaults to no scaling.
  - response_trans_utils`: A named tuple containing to scale/ modify the response.
  - `response_fields` : fields of response to be used for inference
  - `kwargs` : keyword arguments to be splatted into sampling function
"""
function stochastic_inverse(r_obs::resp1, err_resp::resp2, vars, alg_cache::mcmc_cache;
        n_chains=2, model_trans_utils::NamedTuple=(;), # need to take care of this
        response_trans_utils::NamedTuple=(;), params=(;), response_fields=Symbol[],
        kwargs...) where {resp1 <: AbstractResponse, resp2 <: AbstractResponse}
    model_fields = Symbol[]
    # modelD = []
    const_data = []

    # segregate the constants and the Distribution parts of the alg_cache

    apriori = to_dist_nt(alg_cache.apriori)

    for k in keys(apriori)
        if typeof(getfield(apriori, k)) <: Distribution
            push!(model_fields, k)
            push!(const_data, rand(getfield(apriori, k)))
        else
            push!(const_data, getfield(apriori, k))
        end
    end

    likelihood = to_dist_nt(alg_cache.likelihood)
    if response_fields == Symbol[]
        for k in keys(likelihood) # similarly, here it will be propertynames for likelihood being a NamedTuple
            if typeof(getfield(likelihood, k)) <: Function
                push!(response_fields, k)
            end
        end
    end

    # putting trans_utils together for all the fields

    trans_utils_arr = []
    for k in keys(apriori)
        if k in keys(model_trans_utils)
            push!(trans_utils_arr, model_trans_utils[k])
        else
            push!(trans_utils_arr, lin_tf)
        end
    end

    transf_utils = (; zip(keys(apriori), trans_utils_arr)...) # NamedTuple for trans_utils and defaults

    trans_utils_arr = []
    for k in keys(to_resp_nt(r_obs))
        if k in keys(response_trans_utils)
            push!(trans_utils_arr, getfield(response_trans_utils, k))
        else
            push!(trans_utils_arr, lin_tf)
        end
    end

    response_trans_utils = (; zip(keys(to_resp_nt(r_obs)), trans_utils_arr)...)

    @show keys(response_trans_utils)

    m_type = sample_type(alg_cache.apriori)

    if isempty(params)
        params = default_params(m_type)
    end

    # @show to_resp_nt(r_obs)
    # @show to_resp_nt(err_resp)
    msg = """
    variables to be inferred : $(model_fields)
    variables used for inference : $(response_fields)
    model type : $(m_type)
    """
    @info msg

    mcmc_model = mcmc_turing(m_type, const_data, vars, to_resp_nt(r_obs), # ::NamedTuple
        to_resp_nt(err_resp), # ::response
        apriori, # ::NamedTuple
        likelihood, # ::responseDistribution
        params; response_fields=Symbol.(response_fields),
        model_fields=Symbol.(model_fields), model_trans_utils=transf_utils,
        response_trans_utils=response_trans_utils)

    if typeof(alg_cache.sampler).name.module === Pigeons
        n_rounds = Int(round(log2(alg_cache.n_samples)))
        pt = pigeons(; target=TuringLogPotential(mcmc_model), n_chains=n_chains, # Λ ~ 6
            n_rounds=n_rounds,   # low to speed up CI
            record=[traces; round_trip; record_default()], kwargs...)
        return Chains(pt)
    else
        if typeof(alg_cache.sampler) <: Turing.AdvancedVI.VariationalInference
            return vi(mcmc_model, alg_cache.sampler)
        else
            return Turing.sample(mcmc_model, alg_cache.sampler,
                alg_cache.n_samples; verbose=false, kwargs...)
        end
    end
end
