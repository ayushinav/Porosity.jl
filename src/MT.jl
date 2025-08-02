module MT
using InteractiveUtils
using LinearAlgebra
using CairoMakie
using LinearSolve
using NonlinearSolve
using Optimization, OptimizationOptimJL
using Turing
using Distributions
using Pigeons
using Statistics
using SpecialFunctions
using QuadGK
using Enzyme
using DifferentiationInterface
using UnPack
using ProgressMeter
import Base: show

include("abstract_types.jl")
include("probabilistic/init_distributions.jl")
include("inverse/bounds_transformation.jl")
include("geophysics/mt/types.jl")
include("geophysics/mt/forward.jl")
include("geophysics/pretty_printing.jl")

include("rock_physics/conductivity/cache.jl")
include("rock_physics/conductivity/types.jl")
include("rock_physics/conductivity/utils.jl")
include("rock_physics/conductivity/forward.jl")

include("rock_physics/elastic/utils.jl")
include("rock_physics/elastic/cache.jl")
include("rock_physics/elastic/types.jl")
include("rock_physics/elastic/forward.jl")

include("rock_physics/viscous/utils.jl")
include("rock_physics/viscous/cache.jl")
include("rock_physics/viscous/types.jl")
include("rock_physics/viscous/forward.jl")

include("rock_physics/anelastic/cache.jl")
include("rock_physics/anelastic/types.jl")
include("rock_physics/anelastic/utils.jl")
include("rock_physics/anelastic/forward.jl")

include("rock_physics/utils.jl")
include("rock_physics/mixing_phases_core.jl")
include("rock_physics/mixing_phases.jl")
include("rock_physics/combine_models.jl")

include("rock_physics/tune_rp_models.jl")
include("rock_physics/thermals.jl")

include("rock_physics/pretty_printing.jl")

include("probabilistic/rock_physics/conductivity.jl")
include("probabilistic/rock_physics/elasticity.jl")
include("probabilistic/rock_physics/viscosity.jl")
include("probabilistic/rock_physics/anelasticity.jl")
include("probabilistic/rock_physics/mixing_phases.jl")
include("probabilistic/rock_physics/combine_models.jl")
include("probabilistic/rock_physics/tune_rp_models.jl")

include("utils.jl")
include("inverse/utils.jl")
include("inverse/jacobian.jl")
include("inverse/occam.jl")
include("inverse/inv.jl")
include("inverse/nl_inv.jl")
include("inverse/opt_inv.jl")
include("probabilistic/respDistribution.jl")
include("probabilistic/utils.jl")
include("probabilistic/core.jl")
include("probabilistic/inverse.jl")
include("probabilistic/rto.jl")
include("probabilistic/post_inv_utils.jl")
include("plots/utils.jl")
include("plots/plots.jl")
include("plots/prob_utils.jl")

# export μ

# abstracts
export AbstractModel, AbstractResponse
export AbstractGeophyModel, AbstractGeophyResponse
export AbstractModelDistribution, AbstractResponseDistribution
export AbstractGeophyModelDistribution, AbstractGeophyResponseDistribution
export AbstractRockphyModelDistribution, AbstractRockphyResponseDistribution
export AbstractCondModel, AbstractElasticModel, AbstractViscousModel, AbstractAnelasticModel

# geophysics
export MTModel, MTResponse

# rock physics
export RockphyCond, RockphyElastic, RockphyViscous, RockphyAnelastic
export SEO3, UHO2014, Jones2012, Poe2010, Yoshino2009, Wang2006, const_matrix
export Ni2011, Sifre2014, Gaillard2008
export Dai_Karato2009, Zhang2012
export Yang2011
export anharmonic, anharmonic_poro, SLB2005
export HZK2011, HK2003, xfit_premelt
export andrade_psp, eburgers_psp, premelt_anelastic, xfit_mxw, andrade_analytical
export HS1962_plus, HS1962_minus, single_phase, MAL
export HS_plus_multi_phase, HS_minus_multi_phase, GAL
export two_phase_modelType, two_phase_model, multi_phase_modelType, multi_phase_model
export multi_rp_modelType, multi_rp_model, multi_rp_response
export tune_rp_modelType

## thermals
export solidus_Hirschmann2000, solidus_Katz2003
export ΔT_h2o_Katz2003, ΔT_h2o_Blatter2022, ΔT_co2_Dasgupta2007, ΔT_co2_Dasgupta2013
export get_Cco2_m, get_Ch2o_m, get_melt_fraction

# forward
export get_Z, get_appres, get_phase, forward!, forward

#inverse
export transform_utils, default_tf, log_tf, lin_tf, sigmoid_tf, log_sigmoid_tf, pow_tf,
       pow_sigmoid_tf, phi_scale_tf
export mt_jacobian_cache, jacobian_mt, jacobian!
export occam_cache, Occam, nl_cache, NonlinearAlg, opt_cache, OptAlg, linsolve!, occam_step!
export inverse!
export ∂, χ², linear_utils, inverse_utils

# probabilistic
export mcmc_cache, rto_cache
export stochastic_inverse
export normal_dist, uniform_dist

## geophysics
export MTModelDistribution, MTResponseDistribution

## rock physics
export RockphyCondDistribution, RockphyElasticDistribution, RockphyViscousDistribution,
       RockphyAnelasticDistribution
export SEO3Distribution, UHO2014Distribution, Jones2012Distribution, Poe2010Distribution,
       Yoshino2009Distribution, Wang2006Distribution, Dai_Karato2009Distribution,
       Zhang2012Distribution, Yang2011Distribution, const_matrixDistribution
export Ni2011Distribution, Sifre2014Distribution, Gaillard2008Distribution
export anharmonicDistribution, anharmonic_poroDistribution, SLB2005Distribution
export HZK2011Distribution, HK2003Distribution, xfit_premeltDistribution
export andrade_pspDistribution, eburgers_pspDistribution, premelt_anelasticDistribution,
       xfit_mxwDistribution, andrade_analyticalDistribution
export two_phase_modelDistributionType, two_phase_modelDistribution
export multi_phase_modelDistributionType, multi_phase_modelDistribution
export multi_rp_modelDistributionType, multi_rp_modelDistribution,
       multi_rp_responseDistributionType, multi_rp_responseDistribution
export tune_rp_modelDistributionType, tune_rp_modelDistribution

# plots
export plot_response, plot_response!
export plot_model, plot_model!

## post_prob
export get_model_list
export get_kde_image, get_kde_image!, get_mean_std_image, get_mean_std_image!

## utils
export default_params, sample_type

end
