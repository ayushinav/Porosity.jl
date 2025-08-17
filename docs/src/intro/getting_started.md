```@setup rp_intro
using Porosity
```
## Modeling

Modeling for all the functions can be done by calling the `forward` function on the rock physics models.

A very simple example is estimating log conductivity using `SEO3`, a temperature dependent conductivity law for olivine developed by [Constable, 2003](https://doi.org/10.1111/j.1365-246X.2006.03041.x).

Construct the model by passing the needed arguments, here temperature.

```@example rp_intro
T = 1000 + 273.0
model_SEO3 = SEO3(T)
```

and then calling `forward` on it

```@example rp_intro
forward(model_SEO3, [])
```

The `[]` argument is needed to make the interface consistent. It is not used anywhere and you can pass any value to it.

You can pass a whole vector for temperature and get the result :

```@example rp_intro
T = (800:100:1200) .+ 273.0
model_SEO3 = SEO3(T)
forward(model_SEO3, [])
```

## Dimensions

Dimensions for all the rock physics types follow the broadcasting rule, that is, if there are two parameters, one with shape `(11,1)` and other with `(1, 4)` then the output/s will be of `(11, 4)`. This is useful when you want to get the output for multiple values. You can obviously iterate but broadcasting would be usually faster, eg. :

```@example rp_intro
T = (600:100:1600) .+ 273.0
P = 2.0
dg = [1.0, 3.0, 10.0, 30.0]'
σ = 10e-3
ϕ = 0.015

m = HZK2011(T, P, dg, σ, ϕ)
resp = forward(m, []);
println("size of T : ", size(T))
println("size of broadcasted array (same as output): ", size(T .+ P .+ dg .+ σ .+ ϕ)) # 
println("size of ϵ_rate : ", size(resp.ϵ_rate))
println("size of η : ", size(resp.η))
```

Notice that the dimensions of `T` and `dg` were `(11,1)` and `(1, 4)`  respectively and the output shapes are `(11, 4)`.

Ofcourse, simple broadcasting rules dictate that if `T` were of `(11,1)` (or `(11,)`) and so were `dg`, then the output would be `(11, 1)` (or `(11,)`), eg.:

```@example rp_intro
T = (600:100:1600) .+ 273.0
P = 2.0 .+ zero(T) # same size as T
dg = collect(range(0.1; stop=30.0, length=length(T))) # same size as T
σ = 10e-3
ϕ = 0.015

m = HZK2011(T, P, dg, σ, ϕ)
resp = forward(m, []);

println("size of T : ", size(T))
println("size of P : ", size(P))
println("size of broadcasted array (same as output): ", size(T .+ P .+ dg .+ σ .+ ϕ))
println("size of ϵ_rate : ", size(resp.ϵ_rate))
println("size of η : ", size(resp.η))
```

### Anelastic case

For anelastic types, one of the input arguments is frequency `f`. Make sure the  `f` is of the size one more than the dimensions of other values. eg. :

```@example rp_intro
T = (600:100:1400) .+ 273.0
P = 0.2
dg = 3.1
σ = 10.0 * 1.0f-3
ϕ = 0.0
ρ = 3300.0
T_solidus = 1100 + 273.0

f = 10 .^ -collect(range(-2, 4; length=100))'

m = premelt_anelastic(T, P, dg, σ, ϕ, ρ, 0.0, T_solidus, f)
resp = forward(m, [])

println("size of T : ", size(T))
println("size of P : ", size(P))
println(
    "size of broadcasted array without f: ", size(T .+ P .+ dg .+ σ .+ ϕ .+ ρ .+ T_solidus))
println("size of f : ", size(f))
println("size of Qinv (one of the outputs, averaged over frequency) : ", size(resp.Qinv))
```

```@example rp_intro
T = reshape((600:100:1400) .+ 273.0, :, 1, 1)
P = 0.2
dg = reshape(3.0:10.0, 1, :, 1)
σ = 10.0 * 1.0f-3
ϕ = 0.0
ρ = 3300.0
T_solidus = 1100 + 273.0

f = 10 .^ -reshape(collect(range(-2, 4; length=100)), 1, 1, :)

m = premelt_anelastic(T, P, dg, σ, ϕ, ρ, 0.0, T_solidus, f)
resp = forward(m, [])

println("size of T : ", size(T))
println("size of P : ", size(P))
println("size of broadcasted array without f (same as output): ",
    size(T .+ P .+ dg .+ σ .+ ϕ .+ ρ .+ T_solidus))
println("size of f : ", size(f))
println("size of Qinv (one of the outputs) : ", size(resp.Qinv))
println("size of Vave (one of the outputs, averaged over frequency) : ", size(resp.Vave))
```

## Changing initial parameters

A lot of rock physics models have empirical constants on defined physical laws from lab experiments. The `forward` function has a keyword argument `params` which has default value for these constants. These default parameters can be observed by `default_params` function dispatched on the rock physics type, eg:

```@example rp_intro
default_ps_SEO3 = default_params(SEO3)
```

These default values are used by default for all the rock physics models. The `forward` function implicitly calls these default parameters

```@example rp_intro
T = 1000 + 273.0
m = SEO3(T)
forward(m, [])
```

Above is same as :

```@example rp_intro
T = 1000 + 273.0
m = SEO3(T)
forward(m, [], default_ps_SEO3)
```

If you want to use different values for certain laws, you need to create another `params` object, eg. Let's change `S_bfe` and `S_bmg` to different values :

```@example rp_intro
new_params_SEO3 = (; default_ps_SEO3..., S_bfe=1.06f24, S_umg=1.58f26)
```

and then we should see a different output

```@example rp_intro
T = 1000 + 273.0
m = SEO3(T)
forward(m, [], new_params_SEO3)
```

Notice the change in outputs from `forward`.
More complicated examples of changing initial parameters are presented in [Andrade pseudospectral modeling](../models/anelasticity.md#Andrade-pseudoperiod-model) and [Extended Burgers modeling](../models/anelasticity.md#Extended-Burgers-model)

## In-place operations

We do not support in-place operations for rock physics models currently. The rock physics models are generally simple and fast, and the effort to make them allocation-free might not be substantial but watch this space in future for in-place (mutating) functions.
