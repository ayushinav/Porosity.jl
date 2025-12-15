module Porosity
using Reexport
@reexport using SubsurfaceCore
using InteractiveUtils
using LinearAlgebra
using Distributions
using Statistics
using SpecialFunctions
using QuadGK
using UnPack
using NonlinearSolve
import Base: show

import SubsurfaceCore: to_resp_nt, to_nt, forward_helper, forward, sample_type, to_dist_nt,
                       from_nt, default_params

include("models/conductivity/cache.jl")
include("models/conductivity/types.jl")
include("models/conductivity/utils.jl")
include("models/conductivity/forward.jl")

include("models/elastic/utils.jl")
include("models/elastic/cache.jl")
include("models/elastic/types.jl")
include("models/elastic/forward.jl")

include("models/viscous/utils.jl")
include("models/viscous/cache.jl")
include("models/viscous/types.jl")
include("models/viscous/forward.jl")

include("models/anelastic/cache.jl")
include("models/anelastic/types.jl")
include("models/anelastic/utils.jl")
include("models/anelastic/forward.jl")

include("models/utils.jl")
include("models/mixing_phases_core.jl")
include("models/mixing_phases.jl")
include("models/combine_models.jl")

include("models/tune_rp_models.jl")
include("models/thermals.jl")

include("models/pretty_printing.jl")

include("probabilistic/models/conductivity.jl")
include("probabilistic/models/elasticity.jl")
include("probabilistic/models/viscosity.jl")
include("probabilistic/models/anelasticity.jl")
include("probabilistic/models/mixing_phases.jl")
include("probabilistic/models/combine_models.jl")
include("probabilistic/models/tune_rp_models.jl")
include("probabilistic/utils.jl")

include("utils.jl")

include("plots/utils.jl")

# export μ

# abstracts
export AbstractRockphyModelDistribution, AbstractRockphyResponseDistribution
export AbstractCondModel, AbstractElasticModel, AbstractViscousModel, AbstractAnelasticModel

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

## utils
export to_resp_nt, forward_helper, forward, sample_type, to_dist_nt, from_nt
export default_params, sample_type

end
