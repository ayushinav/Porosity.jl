# Model

## Demo

Specifying a model is easy. Let's say you want to add a resistivity distribution with 4 layers as:

  - Layer 1: 500 $\Omega m$, thickness= 100 $m$
  - Layer 2: 1000 $\Omega m$, thickness= 500 $m$
  - Layer 3: 10 $\Omega m$, thickness= 3000 $m$
  - Layer 4: 100 $\Omega m$, half-space

Define the `model` :

```@example model_demo
using MT
ρ = log10.([500.0, 100.0, 400.0, 1000.0]);
h = [100.0, 100.0, 100.0];
m = MTModel(ρ, h)
```

This `model` can then just be passed into `forward` function to get the `response`.

!!! note
    

Always use `Float64` or `Float32` types while defining the vectors for resistivities and thickness. This is done for performance while not imposing any serious constraints since most of the data is generally processed using `Float64` on most CPUs and `Float32` on most GPUs.

The model can then be plotted using

```@example model_demo
f, ax = plot_model(m)
f
```
