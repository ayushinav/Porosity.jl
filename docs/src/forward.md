# Forward modeling

## Demo

Currently, only the recursion solution for MT forward modeling is supported. Once a model is defined, we can get the estimate as:

```@example forward_demo
using MT, CairoMakie
ρ = log10.([500.0, 100.0, 400.0, 1000.0]);
h = [100.0, 1000.0, 10000.0];
m = MTModel(ρ, h)

T = 10 .^ (range(-1, 5; length=57));
ω = 2π ./ T;

resp = forward(m, ω);
nothing # hide
```

## Plots

Since MT `MTResponse` has two fields, and we also want to see the curves for other `MTModel`s, usually for inversion:

```@example forward_demo
f, axs = plot_response(ω, resp; label="1st")
f
```

Another `MTResponse` can be overlain using:

```@example forward_demo
ρ = log10.([100.0, 10.0, 400.0]);
h = [100.0, 10000.0];
m2 = MTModel(ρ, h)
resp2 = forward(m2, ω)
plot_response!(axs, ω, resp2; label="2nd")

f[2, 2] = Legend(f, axs[1])
f
```

Note that we've plotted the second response curve as a line, which can simply be achieved by passing `plt_type = :plot`

## Benchmark

In-place operations for fast non-allocating computations are also supported. These would greatly speed up the inverse processes.

```@example forward_demo
using BenchmarkTools
@benchmark forward!(resp, m, ω)
```

The above benchmark was done on Mac M1.

If we want, we can get the data on a different scale ((using domain transformation)[domain_transformation.md]) than is provided by default using `:

```@example forward_demo
ρ = log10.([500.0, 100.0, 400.0, 1000.0]);
h = [100.0, 1000.0, 10000.0];
m = MTModel(ρ, h)

T = 10 .^ (range(-1, 5; length=57));
ω = 2π ./ T;

resp = forward(m, ω, (ρₐ=log_tf, ϕ=lin_tf));
f, ax = plot_response(ω, resp; label=false)
f
```
