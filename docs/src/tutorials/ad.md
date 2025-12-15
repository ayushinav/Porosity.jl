# Automatic Differentiation

```@setup example_ad
using Porosity
```

## Example

We can compute derivatives of the models used in the package using AD ([automatic differentiation](https://en.wikipedia.org/wiki/Automatic_differentiation)). While the functionalities that require AD are built-in, we demonstrate below how the derivatives can be explicitly obtained. We take the example of [`UHO2014`](@ref) which estimates conductivity of olivine depending on temperature and water content.

Let's first define the model

```@example example_ad
m = UHO2014(1000.0, 100.0)
```

```@example example_ad
resp = forward(m, [])
```

For AD, we use the [`DifferentiationInterface`](https://juliadiff.org/DifferentiationInterface.jl/DifferentiationInterface/stable/) library with [`Enzyme`](https://enzyme.mit.edu/julia/stable/) backend.

```
using DifferentiationInterface, Enzyme

f(T) = forward(UHO2014(T, 100.) , []).σ
ad_wrt_T = derivative(f, AutoEnzyme(), 1000.)
```

Let's verify that using finite differences :

```
fd_wrt_T = 0.5 * (f(1001.) - f(999.))
```

Similarly, we can do the same for water content :

```
f(Cw) = forward(UHO2014(1000., Cw) , []).σ
ad_wrt_Cw = derivative(f, AutoEnzyme(), 100.)
```

```
fd_wrt_Cw = 0.5 * (f(101.) - f(99.))
```

If you can imagine the input as a vector, you can obtain the  gradient to get the derivative in one step :

```
function f(p)
    T = p[1]
    Cw = p[2]
    forward(UHO2014(T, Cw) , []).σ
end

ad_p = DifferentiationInterface.gradient(f, AutoEnzyme(), [1000., 100.])
```

::: info
If you need derivatives explicitly, it is recommended to check out `DifferentiationInterface.jl` more thoroughly, especially the mutating counterparts (`derivative!` and `gradient!`).
:::

## A bit of caution

Please note that you would not need to compute the derivatives in most cases. If you are using a neural network you can just plug the forward model in front of the network to train it. This page is to show the compatibility of the code for automatic differentiation.

An MWE (minimal working example) for such a model would be :

```julia
using Porosity, Lux

function f(p)
    T = p[1, :]
    Cw = p[2, :]
    forward(UHO2014(T, Cw), []).σ
end

nn = Chain(Dense(1 => 4, tanh), Dense(4 => 4, tanh), Dense(4 => 2), f)
```

The above network can then be trained by using the same instructions as for training a `Lux` network.
