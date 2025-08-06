mutable struct PseudoGeophysicalModel{
    T1 <: AbstractArray{<:Any}, T2 <: AbstractArray{<:Any}} <: AbstractGeophyModel{T1, T2}
    m::T1
    h::T2
end
