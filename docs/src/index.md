```@raw html
---
# https://vitepress.dev/reference/default-theme-home-page
layout: home

hero:
  name: Porosity.jl
  text: Rock physics parameter inference in Julia
  tagline: Automatic Differentiation enabled probabalistic inference of rock physics
  image:
    src: logo.png
    alt: Porosity.jl
  actions:
    - theme: alt
      text: View on Github
      link: https://github.com/ayushinav/Porosity.jl
    
features:
  - icon: 🔢
    title: Modeling
    details: Estimate geophysical observables using rock physics
    link: intro/getting_started

  - icon: 📊
    title: Probabilstic inference
    details: Perform probablistic inference of parameters
    link: /tutorials/stochastic_inverse

  - icon: ∂
    title: Differentiability
    details: Get derivatives using automatic differntiation
    link: /tutorials/ad
---
```

## Installation

You can install `Porosity.jl` on Julia by running:

```julia
using Pkg
Pkg.add("Porosity.jl")
```

A large part of the package has functionalities borrowed from the [Very Broadband Rheology Calculator (VBRc)](https://github.com/vbr-calc/vbr), coded up in MATLAB.

!!! tip "Very Broadband Rheology Calculator"

    MATLAB framework for calculating material properties from thermodynamic state variables (e.g., temperature, pressure, melt fraction, grain size) using a wide range of experimental scalings.


## Available models

::: code-group

``` [Conductivity]
* Dai_Karato2009
* Gaillard2008
* Jones2012
* Ni2011
* Poe2010
* SEO3
* Sifre2014
* UHO2014
* Wang2006
* Yang2011
* Yoshino2009
* Zhang2012
* const_matrix

The whole list can be obtained by running `subtypes(AbstractCondModel)`
```

``` [Elastic]
* SLB2005
* anharmonic
* anharmonic_poro

The whole list can be obtained by running `subtypes(AbstractElasticModel)`
```

``` [Viscous]
* HK2003
* HZK2011
* xfit_premelt

The whole list can be obtained by running `subtypes(AbstractViscousModel)`
```

``` [Anelastic]
* andrade_analytical
* andrade_psp
* eburgers_psp
* premelt_anelastic
* xfit_mxw

The whole list can be obtained by running `subtypes(AbstractAnelasticModel)`
```

:::