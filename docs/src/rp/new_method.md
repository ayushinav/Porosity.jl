# Contributing

## Adding a new rock physics equation

If you want to add a new rock physics equation, eg., a conductivity mechanism that got out recently and fits with your research or an anelasticity model you want to add, please open an issue/ pull request.

You can contribute by adding a new method following the steps listed below where we follow how a `conductivity` model is added but similar idea is used for `viscous`/`elastic`/`anelastic`:

 1. add a new `mutable struct` similar to other fomulations in `src/conductivity/types.jl`
 2. define `forward` dispatch in `src/conductivity/forward.jl`
 3. add the constants to be used in `forward` dispatch to `src/conductivity/cache.jl`. It is recommended to look at other `params` defined in the `cache.jl` file. Adding empirical constants in `params` keeps the codebase cleaner and organized.
 4. Define dispatch on `default_params` in `src/conductivity/types.jl` itself.
 5. If required, add other functions in `src/conductivity/utils.jl`
 6. add a new distribution struct in `probabilistic/conductivity.jl` similar in definition as the `struct` defined in step 1.
 7. Finally export from the `src/MT.jl` file.

It's always a good idea to add some docstring and test for the methods you add.

## Adding a new phase mixing method

To add a new phase mixing method, we follow the example of [`HS_1962_plus`](@ref) :

 1. Define a new `mutable struct` in `src/mixing_phases.jl`
 2. Define the function dispatch `mix_model` in `src/mixing_phases_core.jl` which dictates how bulk property can be estimated.
 3. Define dispatch on `from_nt` which basically outputs the object of the same type. This was needed to work things out for `MAL` and `GAL`.
 4. Finally export from the `src/MT.jl` file.

Please note that the 2-phase mixing models are different from a more general N-phase mixing model. There isn't much difference when adding a new method except that if it's a multi-phase mixing model, add that to `multi_phase_mix_types`, which is a union of all types which work on [`multi_phase_model`](@ref) whereas if it's a two-phase model, add that to `two_phase_mix_types`.
