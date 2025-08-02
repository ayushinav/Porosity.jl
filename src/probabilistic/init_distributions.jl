"""
    mutable struct modelDistribution{T1<: Union{Distribution, AbstractArray}, T2<: Union{Distribution, AbstractArray}} # where T1,T2 
        m::T1
        h::T2
    end

create a placeholder to store the `Distributions.jl` sampler for a priori
"""
mutable struct MTModelDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray}} <:
               AbstractGeophyModelDistribution
    m::T1
    h::T2
end

"""
    struct responseDistribution{T1<: Union{Function, Nothing}, T2<: Union{Function, Nothing}} # where T1,T2
        ρₐ::T1
        ϕ::T2
    end

create a placeholder to store functions to obtain `Distributions.jl` samplers for the likelihood function
"""
struct MTResponseDistribution{
    T1 <: Union{Function, Nothing}, T2 <: Union{Function, Nothing}} <:
       AbstractGeophyResponseDistribution
    ρₐ::T1
    ϕ::T2
end

"""
    RockphyModelDistribution{
        T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray}} <:
                AbstractRockphyModelDistribution
        params::T1 # vector of parameters 
        p_names::Vector{<:Symbol} # Vector of symbols telling the parameters in vector 
        ϕ::T2 # phase ratios
        model_list::Vector{<:Type}
        mixing_type::Vector

Initializes a placeholder for all the variables to be initialized for use in the `stochastic_inverse`(@ref stochastic_inverse).
Has the same structure as `construct_mixing_models`(@ref construct_mixing_models) but instead inputs distributions, and therefore respects the same form.

## Arguments

  - `params` : A distribution from `Distributions.jl` to sample the rock physics parameters.
  - `p_names` : The `params` vector contains the distributions for the parameters specified by the vector of symbols called `p_names`
  - `ϕ` : Vol. fraction of different phases (can be a distribution or a vector)
  - `model_list` : list of models to be used in mixing
  - `mixing_type` : Type of mixing model to be used
"""
# mutable struct RockphyModelDistribution{
#     T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray}} <:
#                AbstractRockphyModelDistribution
#     params::T1 # vector of parameters 
#     p_names::Vector{<:Symbol} # Vector of symbols telling the parameters in vector 
#     ϕ::T2 # phase ratios
#     model_list::Vector{<:Type}
#     mixing_type::Vector
#     # add water and partition ratios

#     # @assert ϕ and model_lists have same length
#     # @assert p contains all the variables required by all models in model_list
# end

# struct RockphyResponseDistribution{T} <: AbstractResponseDistribution
#     σ::T
# end
