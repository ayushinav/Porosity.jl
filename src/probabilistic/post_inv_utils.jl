"""
    get_model_list(chains::chain, mDist::mdist; 
        trans_utils = (m = pow_tf, h = lin_tf,)) where {mdist <: AbstractModelDistribution, chain <: Chains}

returns a list of models from the `Chains` variable obtained from [`stochastic_inverse`](@ref)

## Arguments

  - `chains` : `Chains` object obtained from the `Turing` model
  - `mDist` : *a priori* distribution defined before performing stochastic inversion
"""
function get_model_list(chains::chain,
        mDist::mdist;
        trans_utils=(m=lin_tf, h=lin_tf)) where {
        mdist <: AbstractModelDistribution, chain <: Chains}
    model_type = sample_type(mDist)

    length_vec = zeros(Int64, length(propertynames(mDist)))
    for (i, k) in enumerate(propertynames(mDist))
        length_vec[i] = (typeof(getproperty(mDist, k)) <: Distribution) ?
                        length(rand(getproperty(mDist, k))) : length(getproperty(mDist, k))
    end

    preds = []
    for k in chains.name_map.parameters
        push!(preds, chains[k].data[:])
    end
    pred = hcat(preds...)
    model_list = []

    for i in 1:size(pred, 1)
        m = []
        prev_length = 1
        for (j, k) in enumerate(fieldnames(model_type))
            if typeof(getproperty(mDist, k)) <: Distribution
                next_length = prev_length + length_vec[j]
                push!(m, pred[i, prev_length:(next_length - 1)])
                prev_length = 0 + next_length
            else
                push!(m, getproperty(mDist, k))
            end
        end

        model_sample = model_type(m...)

        for k in fieldnames(model_type)
            if k in keys(trans_utils)
                getfield(model_sample, k) .= broadcast(
                    getproperty(trans_utils[k], :tf), getfield(model_sample, k))
            end
        end

        push!(model_list, model_sample)
    end
    return model_list
end

"""
    get_ρ_at_z(pred, zs)

returns the values of the model splatted as a vector in `pred` at points defined by `zs`
"""
function get_ρ_at_z(pred, zs)
    n = 1 + length(pred) ÷ 2
    h = cumsum(pred[(n + 1):end])
    res = zeros(length(zs))
    idx = zs .<= h[1]
    res[idx] .= pred[1]

    for i in 2:length(h)
        idx = ((h[i - 1] .< zs .<= h[i]))
        res[idx] .= pred[i]
    end
    idx = zs .> h[end]
    res[idx] .= pred[n]
    return res
end
