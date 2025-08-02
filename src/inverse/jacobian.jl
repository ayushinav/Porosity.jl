# """
# `jacobian_cache`: cache `struct` to store the values during the estimation of jacobian.
# """
# mutable struct jacobian_cache{T1, T2, T3}
#     r1::T1
#     r2::T2
#     j::T3
# end

# function jacobian_cache(resp_fields, resp, mod, mod_fields)
#     nresps = sum([length(getfield(resp, k)) for k in resp_fields])
#     nmods = sum([length(getfield(mod, k)) for k in mod_fields])

#     eltype_ = eltype(getfield(mod, first(mod_fields)))

#     return jacobian_cache(zero(resp), zero(resp), zeros(eltype_, nresps, nmods))#resp_type(arr...))
# end

# """
#     jacobian!(m::model, 
#         vars::Vector{T},
#         jc::jacobian_cache;
#         model_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(m))], 
#         response_fields::Vector{Symbol}= [k for k ∈ fieldnames(typeof(J))], 
#         response_trans_utils::NamedTuple=(; ρₐ=lin_tf, ϕ=lin_tf)
#         ) where T <: Union{Float64, Float32}

# overwrites a `jacobian_cache` cache to calculate the jacobian of a `model`.
# If no `model_fields` or `response_fields` are passed, all the fields of `model` and the `response` (defined in `jacobian`) will be used.
# """
# function jacobian!(jc::jacobian_cache,
#         m::model,
#         vars::Vector{T},
#         model_fields::Vector{Symbol},
#         response_fields::Vector{Symbol},
#         response_trans_utils::resp_trans_utils_T=MT.default_mt_tf_fns) where {
#         T <: Union{Float64, Float32}, model <: AbstractModel, resp_trans_utils_T}
#     n_mod_start = 0
#     ϵ = sqrt(eps(eltype(getfield(m, first(model_fields)))))
#     @inbounds for k in model_fields
#         ϵ = sqrt(eps(eltype(getfield(m, k))))

#         @inbounds for i in eachindex(getfield(m, k))
#             getfield(m, k)[i] = getfield(m, k)[i] + ϵ
#             forward!(jc.r1, m, vars, response_trans_utils)

#             getfield(m, k)[i] = getfield(m, k)[i] - 2ϵ
#             forward!(jc.r2, m, vars, response_trans_utils)

#             n_resp_start = 1
#             n_resp_end = 0

#             @inbounds for l in response_fields
#                 n_resp_end += length(getfield(jc.r1, l))
#                 view(jc.j, n_resp_start:n_resp_end, n_mod_start + i) .= getfield(
#                     jc.r1, l) .- getfield(jc.r2, l)

#                 n_resp_start += n_resp_end
#             end
#             getfield(m, k)[i] = getfield(m, k)[i] + ϵ
#         end
#         n_mod_start += length(getfield(m, k))
#     end
#     rmul!(jc.j, inv(2ϵ))
#     nothing
# end

# # jacobian works fine for 1D MT model, we would have to check for 2D/3D models.
