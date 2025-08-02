
# models

"""
Abstract model type that is the supertype of all `model`s in the package.
"""
abstract type AbstractModel end

"""
Abstract model type that is the supertype of all geophysical models in the package.
"""
abstract type AbstractGeophyModel <: AbstractModel end

"""
Abstract model type that is the supertype of all rock models in the package.
"""
abstract type AbstractRockphyModel <: AbstractModel end

# responses

"""
Abstract model type that is the supertype of all `response`s in the package.
"""
abstract type AbstractResponse end

"""
Abstract model type that is the supertype of all geophysical responses in the package.
"""
abstract type AbstractGeophyResponse <: AbstractResponse end

"""
Abstract model type that is the supertype of all rock physics responses in the package.
"""
abstract type AbstractRockphyResponse <: AbstractResponse end

# distributions

## modes

"""
Abstract model type that is the supertype of all `model` distributions in the package.
"""
abstract type AbstractModelDistribution end

"""
Abstract model type that is the supertype of all geophysical `model` distributions in the package.
"""
abstract type AbstractGeophyModelDistribution <: AbstractModelDistribution end

"""
Abstract model type that is the supertype of all rock physics `model` distributions in the package.
"""
abstract type AbstractRockphyModelDistribution <: AbstractModelDistribution end

## responses

"""
Abstract model type that is the supertype of all `response` distributions in the package.
"""
abstract type AbstractResponseDistribution end

"""
Abstract model type that is the supertype of all geophysical `response` distributions in the package.
"""
abstract type AbstractGeophyResponseDistribution <: AbstractResponseDistribution end

"""
Abstract model type that is the supertype of all rock physics `response` distributions in the package.
"""
abstract type AbstractRockphyResponseDistribution <: AbstractResponseDistribution end
