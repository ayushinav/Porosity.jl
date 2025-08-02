# `MT.jl`

`MT.jl` is supposed to be a high performance code for doing forward and inverse modeling in geophysics using julia. We hope to write the code structure such that any other geophysical survey can also be used and we can tend towards a joint forward and inverse modeling library.

## Forward modeling

While forward modeling typically requires solving a PDE obtained using the quasi-static approximation, in 1D, we are fortunate to have the solution for surface impedance in a more analytical form. Currently, this is what is supported.

Supported methods:

  - 1D Magnetotellurics (MT)

## Inverse modeling

No surprises here that we are almost always trying to solve for an under-determined system.

Deterministic schemes supported:

  - Occam
  - Nonlinear schemes using NonlinearSolve.jl
  - Nonlinear schemes using Optimization.jl

Probabilistic schemes supported:

  - MCMC with fixed grids
  - MCMC with flexible grids
  - RTO-TKO

## Rock physics

We support multiple rock physics models :

  - Conductivity models
  - Elasticity models
  - Viscosity models
  - Anelasticity models
