# COV_EXCL_START
mutable struct PseudoGeophysicalModel{
    T1 <: AbstractArray{<:Any}, T2 <: AbstractArray{<:Any}} <: AbstractGeophyModel{T1, T2}
    m::T1
    h::T2
end

mutable struct PseudoGeophysicalModelDistribution{
    T1 <: Union{Distribution, AbstractArray}, T2 <: Union{Distribution, AbstractArray}} <:
               AbstractGeophyModelDistribution{T1, T2}
    m::T1
    h::T2
end

SubsurfaceCore.sample_type(d::PseudoGeophysicalModelDistribution) = PseudoGeophysicalModel
# COV_EXCL_STOP
