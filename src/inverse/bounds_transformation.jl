"""
`sigmoid(m)`:
move to model domain from optimization domain using a sigmoid transformation
"""

function sigmoid(m::T1, bounds::T2) where {T1, T2} #m::T1, bounds::T2) where {T1, T2} #m::T1, bounds::Vector{T})::T where {T <: Union{Float32, Float64}, T1}
    σ(x) = 1 / (1 + exp(-x / (bounds[2] - bounds[1])))
    # return 10^(σ(m)*log10(bounds[2]/bounds[1])+ log10(bounds[1]))
    return σ(m) * (bounds[2] - bounds[1]) + bounds[1]
end

function pow_sigmoid(m::T1, bounds::T2) where {T1, T2} #m::T1, bounds::Vector{T})::T where {T <: Union{Float32, Float64}, T1}
    return exp10(sigmoid(m, bounds))
end

scale_fn(m::T1, scale::T2) where {T1, T2} = m / scale
"""
`d_sigmoid(m)`:
gradient for the transformation from optimization domain to model domain. Used for estimating jacobians, but is also useful in analysing sensitivities.
"""
function d_sigmoid(m::T1, bounds::T2) where {T1, T2} #m::T1, bounds::Vector{T})::T where {T <: Union{Float32, Float64}, T1}
    d_σ(x) = inv((1 + exp(x)) * (1 + exp(-x)))
    # return σ(log10(m))*(log10(bounds[2])- log10(bounds[1]))+ log10(bounds[1])
    return d_σ(m / (bounds[2] - bounds[1]))#*(bounds[2]- bounds[1])
end

function d_pow_sigmoid(m::T1, bounds::T2) where {T1, T2} #m::T1, bounds::Vector{T})::T where {T <: Union{Float32, Float64}, T1}
    return exp10(sigmoid(m, bounds)) * d_sigmoid(m, bounds)
end

d_scale_fn(m::T1, scale::T2) where {T1, T2} = inv(scale)

"""
`inverse_sigmoid()`: get back to the optimization domain from model domain
"""
function inverse_sigmoid(
        x::T1, bounds::Vector{T})::T where {T <: Union{Float32, Float64}, T1}
    # if x ≈ bounds[2] return 100 end
    # if x ≈ bounds[1] return -100 end

    return (bounds[2] - bounds[1]) * (log(abs(x - bounds[1])) - log(abs(bounds[2] - x)))
end

function inverse_pow_sigmoid(m::T1, bounds::T2) where {T1, T2} #m::T1, bounds::Vector{T})::T where {T <: Union{Float32, Float64}, T1}
    return inverse_sigmoid(log10(m), bounds)
end

inverse_scale_fn(x, scale) = x * scale

struct transform_utils{T}
    p::Vector{T} # parameters of the transformation function
    tf::Function
    itf::Function
    dtf::Function
end

"""
    `transform_utils(p, tf, itf, dtf)`

Contains the parameters and functions for transformation from optimization to model domains.

## Arguments

  - `p` : Vector for parameterization, eg. upper and lower bounds of sigmoid transformation
  - `tf` : function to convert from computational to model domain
  - `itf` : function to convert from model to computational domain
  - `dtf` : gradient of `tf` wrt to computational domain variable

## Usage

Implementation of `pow_tf` : when computational domain = log10 (model domain)

```julia
transform_utils(
    Vector{Float32}([]), (x, p) -> 10^x, (x, p) -> log10(x), (x, p) -> (10^x * log(10)))
```

Also checkout relevant documentation
"""
function transform_utils(p::Vector{T}, tf::Function, itf::Function,
        dtf::Function) where {T <: Union{Float32, Float64}}
    return transform_utils{eltype(p)}(
        p, (x) -> tf(x, p), (x) -> itf(x, p), (x) -> dtf(x, p))
end

# should generally be good for most inversions

"""
    sigmoid_tf

A [`transform_utils`](@ref) constant using `σ(x)` ; `σ(x) = 1 / (1 + exp(-x))` \n

  - computational domain to model domain : `m = -3 + 9σ(x)` \n
  - model domain to computational domain : `x = σ⁻¹((m+3)/9)` \n

The above bounds the model domain in [-3, 6], to bound the model parameters in a different domain [a,b] : \n
\n
`transform_utils([a, b], sigmoid, inverse_sigmoid, d_sigmoid)`
"""
const sigmoid_tf = transform_utils([-3.0, 6.0], sigmoid, inverse_sigmoid, d_sigmoid);

"""
    pow_tf

A [`transform_utils`](@ref) constant using `exp10` \n

  - computational domain to model domain : m = exp10(x) \n
  - model domain to computational domain : x = log10(m) \n
"""
const pow_tf = transform_utils(
    Vector{Float32}([]), (x, p) -> exp10(x), (x, p) -> log10(x), (x, p) -> (10^x * log(10)));

"""
    log_tf

A [`transform_utils`](@ref) constant using `log10` \n

  - computational domain to model domain : m = log10(x) \n
  - model domain to computational domain : x = exp10log10(m) \n
"""
const log_tf = transform_utils(
    Vector{Float32}([]), (x, p) -> log10(x), (x, p) -> exp10(x), (x, p) -> inv(x * log(10)));

"""
    pow_sigmoid_tf

A [`transform_utils`](@ref) constant similar to `sigmoid_tf`(@ref) but further performs `exp10` to convert into model domain \n

  - computational domain to model domain : `m = exp10(-3 + 9σ(x))` \n
  - model domain to computational domain : x` = σ⁻¹((log10(m)+3)/9)` \n

The above bounds the model domain in [10⁻³, 10⁶], to bound the model parameters in a different domain [10ᵃ,10ᵇ] : \n
\n
`transform_utils([a, b], sigmoid, inverse_sigmoid, d_sigmoid)`
"""
const pow_sigmoid_tf = transform_utils(
    [-3.0, 6.0], pow_sigmoid, inverse_pow_sigmoid, d_pow_sigmoid);

"""
    lin_tf

A [`transform_utils`](@ref) constant that doesn't transform \n
Default `transform_utils` if nothing else is provided

  - computational domain to model domain : m = x \n
  - model domain to computational domain : x = m \n
"""
const lin_tf = transform_utils(Vector{Float32}([]), (x, p) -> x, (x, p) -> x, (x, p) -> 1.0);

"""
    phi_scale_tf

A [`transform_utils`](@ref) that normalizes linearly on a scale of 90. \n

  - computational domain to model domain : m = x/90 \n
  - model domain to computational domain : x = 90m \n
    The above scales wrt 90, to scale with a different value `a` : \n
    \n
    `phi_scale_tf = transform_utils([a], (x, p) -> scale_fn(x, first(p)), (x, p) -> inverse_scale_fn(x, first(p)), (x, p) -> d_scale_fn(x, first(p)))`
"""
const phi_scale_tf = transform_utils([90.0f0], (x, p) -> scale_fn(x, first(p)),
    (x, p) -> inverse_scale_fn(x, first(p)), (x, p) -> d_scale_fn(x, first(p)))
