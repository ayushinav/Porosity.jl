mutable struct MTResponse{T1, T2} <: AbstractGeophyResponse
    ρₐ::T1
    ϕ::T2
end

#=
More often than not, T1 and T2 will be same, that is typeof(ρₐ) and typeof(ϕ) will be same but this might become different when using a rock physics model in front of this. 
For now, we go with the same design as for `MTModel`.
=#

"""
create a model type for a given resistivity distribution
that can be used to calculate forward response for 1d MT

## Usage

```jldoctest
m = [4.0, 2.0, 3.0]
h = [1000.0, 100.0]
model = MTModel(m, h)
print(model)

# output

1D MTModel : 
Layer \t log(ρ)  h
____________________________
1 \t 4.0 \t 1000.0
2 \t 2.0 \t 100.0
3 \t 3.0 \t ∞
```
"""
mutable struct MTModel{T1 <: AbstractArray{<:Any}, T2 <: AbstractArray{<:Any}} <:
               AbstractGeophyModel
    m::T1
    h::T2
end

# explain in blog the reasoning behind this! This covers all 1D, 2D, 3D models for MT.
# make pretty tables to print these models

default_params(::Type{T}) where {T <: MT.AbstractGeophyModel} = (;)
