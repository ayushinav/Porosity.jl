# Visualization

We make use of `Makie.jl` to generate images. This allows us to be more flexible with and also generate more fancy figures.

Here, we want to give a slightly detailed outlook on generating figures. `Makie.jl` provides multiple backends to generate graphics. We use `CairoMakie` because it uses CPUs and is compatible with most machines. There is also `GLMakie` which makes use of GPUs, and the same functions should work with any backend.

## Using `Makie`

`Makie.jl` works around three objects : `scene`, `figure` and `axis`. Of these three, `scene` is mostly useful for 3D rendering, which, while is a nice feature, is not really useful for us. We, therefore, make use of the latter two.

While a more interested reader is encouraged to explore the `Makie.jl` documentation, here we suffice by saying that a `figure` is the whole figure, containing multiple `axis` objects. Each `axis` can contain a plot, colorbar, legend, etc. Therefore, we advised that while creating legends, use another `axis`. This allows more control.

For the purpose of this tutorial, know that when you call a plotting function, it returns a `figure` and all the `axis` in the figure. You must have guessed by now, but we want to make it slightly more explicit that most of the plotting happens on `axis` while `figure` is just a container around it to control the aspects (pun intended) of the figure. Therefore, when overlaying multiple plots on the same graph, we only need `axis` and not necessarily `figure`. All the mutating functions, therefore, ask for `axis`.

If you want to make an empty figure, just call:

```@example viz_tut
using CairoMakie
f = Figure(; size=(200, 200))
f
```

We now have the container. Time to add an `axis`:

```@example viz_tut
ax = Axis(f[1, 1])
f
```

Let's add another

```@example viz_tut
ax1 = Axis(f[2, 1])
f
```

and another...

```@example viz_tut
ax2 = Axis(f[2, 2])
f
```

Notice the arrangement of axes according to the indices passed for `f` in `Axis`. This allows us to create more fanciful layouts:

```@example viz_tut
f = Figure(; size=(200, 200))
ax1 = Axis(f[1:2, 1])
ax2 = Axis(f[1:2, 2])
ax3 = Axis(f[3, 1])
ax4 = Axis(f[3, 2])
ax5 = Axis(f[1:3, 3])
f
```

The usage will become more apparent with the following example:

## Plotting models

Let's say you have a model defined as :

```@example viz_tut
using MT

ρ = log10.([500.0, 100.0, 400.0, 1000.0]);
h = [100.0, 1000.0, 10000.0];
m = MTModel(ρ, h)
```

Now, we plot the model with :

```@example viz_tut
f, ax = plot_model(m)
f
```

See how the `plot_model` function returns a `figure` and an `axis`. If we want to overlay another model on the same plot, we make use of `plot_model!`

```@example viz_tut
ρ = log10.([50.0, 100.0, 1000.0, 200.0]);
h = [1000.0, 1000.0, 10000.0];
m2 = MTModel(ρ, h)

plot_model!(ax, m)
```

We only needed to pass `ax` to plot on the same `axis`. This might also seem a bit intuitive as well. At this point, we want to encourage creating empty `figure` beforehand and then always using `plot_model!(...)`. This gives us more control on how the axes are arranged among other stuff. Here's an example to make things clear.

```@example viz_tut

ρ = log10.([500.0, 100.0, 400.0, 1000.0]);
h = [100.0, 1000.0, 10000.0];
m = MTModel(ρ, h)

f = Figure(; size=(200, 200)) # creating empty figure
ax = Axis(f[1, 1])

ρ = log10.([50.0, 100.0, 1000.0, 200.0]);
h = [1000.0, 1000.0, 10000.0];
m2 = MTModel(ρ, h)

plot_model!(ax, m)
plot_model!(ax, m2)
f
```

Notice how the above plots did not have legends. Remember that legends were another `axis`.

```@example viz_tut
ρ = log10.([500.0, 100.0, 400.0, 1000.0]);
h = [100.0, 1000.0, 10000.0];
m = MTModel(ρ, h)

f = Figure(; size=(200, 200)) # creating empty figure
ax = Axis(f[1, 1])

ρ = log10.([50.0, 100.0, 1000.0, 200.0]);
h = [1000.0, 1000.0, 10000.0];
m2 = MTModel(ρ, h)

plot_model!(ax, m; label="m1")
plot_model!(ax, m2; label="m2")

Legend(f[1, 2], ax)
f
```

## Plotting responses

Plotting responses is pretty similar.

```@example viz_tut
T = 10 .^ (range(-1, 5; length=57));
ω = 2π ./ T;

resp = forward(m, ω);
f, axs = plot_response(ω, resp; label="1st", plt_type=:scatter)
f
```

`axs` now has two axes because `MTResponse` has two responses, apparent resistivity and phase.
Overlaying is also nothing different. Pass the `axis` in the mutating form of function.

```@example viz_tut
resp2 = forward(m2, ω);
plot_response!(axs, ω, resp; label="2nd", plt_type=:scatter)
f
```

We again want to emphasize that creating an empty figure first followed by the usage of mutating functions is encouraged and gives a lot more control. We again demonstrate the same by including legend in the plots.

```@example viz_tut
f = Figure()
ax1 = Axis(f[1, 1])
ax2 = Axis(f[2, 1])

axs = [ax1, ax2]
plot_response!(axs, ω, resp; label="1st", plt_type=:scatter)
plot_response!(axs, ω, resp2; label="2nd", plt_type=:scatter)

Legend(f[1:2, 2], ax2)
f
```

Creating a `figure` beforehand and following the recommended guidelines also allows us to plot figures as we like:

```@example viz_tut
f = Figure()
ax1 = Axis(f[1:3, 1]; xminorticksvisible=true, backgroundcolor=(:gray, 0.1))
ax2 = Axis(f[1:3, 2]; xminorticksvisible=true, backgroundcolor=(:magenta, 0.1))

axs = [ax1, ax2]
plot_response!(axs, ω, resp; label="1st", plt_type=:scatter)
plot_response!(axs, ω, resp2; label="2nd", plt_type=:scatter)

Legend(f[4, 1:2], ax2; orientation=:horizontal)
f
```

## Concluding remarks

This tutorial is to give you a walkthrough of how `Makie` works in our context. Most other plotting functions in the package follow the same style. There are certain other keyword arguments specific to each of these functions, eg. `half_space_thickness` in `plot_model!()` and `plt_type` in `plot_response!()`. Information on these arguments specific to such functions are available in their respective docstrings.

`Makie.jl` also has a number of these keyword arguments, which can be used to manipulate the figure. Do checkout the respective documentations.
