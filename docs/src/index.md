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

A large part of the package has functionalities borrowed from the [Very Broadband Rheology Calculator (VBRc)](https://vbr-calc.github.io/vbr/), coded up in MATLAB.

```@raw html
<!--
  GitHub-Style Repository Box
  Copy and paste this entire block into your Documenter.jl @raw html block.
  You can customize the link, repository name, and description below.
-->
<div style="padding: 16px; border: 1px solid #d0d7de; border-radius: 6px; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Helvetica, Arial, sans-serif, 'Apple Color Emoji', 'Segoe UI Emoji'; max-width: 80%; margin: 16px 0; position: relative;">
  
  <!-- Icons Container -->
  <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 12px; padding: 0 4px;">
    <!-- GitHub Icon (Right) -->
    <a href="https://github.com/vbr-calc/vbr" target="_blank" rel="noopener noreferrer" style="display: inline-block; text-decoration: none; color: inherit;">
      <svg height="24" width="24" viewBox="0 0 16 16" version="1.1" aria-hidden="true">
        <path fill-rule="evenodd" fill="currentColor" d="M8 0C3.58 0 0 3.58 0 8c0 3.54 2.29 6.53 5.47 7.59.4.07.55-.17.55-.38 0-.19-.01-.82-.01-1.49-2.01.37-2.53-.49-2.69-.94-.09-.23-.48-.94-.82-1.13-.28-.15-.68-.52-.01-.53.63-.01 1.08.58 1.23.82.72 1.22 1.87.87 2.33.66.07-.52.28-.87.51-1.07-1.78-.2-3.64-.89-3.64-3.95 0-.87.31-1.59.82-2.15-.08-.2-.36-1.02.08-2.12 0 0 .67-.21 2.2.82.64-.18 1.32-.27 2-.27.68 0 1.36.09 2 .27 1.53-1.04 2.2-.82 2.2-.82.44 1.1.16 1.92.08 2.12.51.56.82 1.27.82 2.15 0 3.07-1.87 3.75-3.65 3.95.29.25.54.73.54 1.48 0 1.07-.01 1.93-.01 2.2 0 .21.15.46.55.38A8.013 8.013 0 0016 8c0-4.42-3.58-8-8-8z"></path>
      </svg>
    </a>
  </div>

  <div style="display: flex; align-items: center; margin-bottom: 8px;">
    <!-- Repository Name and Link -->
    <a href="https://github.com/vbr-calc/vbr" target="_blank" rel="noopener noreferrer" style="font-weight: 600; font-size: 16px; color: #0969da; text-decoration: none;">
      vbr-calc/vbr
    </a>
  </div>
  <!-- Repository Description -->
  <p style="font-size: 14px; color: #57606a; margin: 0;">
    The Very Broadband Rheology Calculator (VBRc) provides a useful framework for calculating material properties from thermodynamic state variables (e.g., temperature, pressure, melt fraction, grain size) using a wide range of experimental scalings. The VBRc at present contains constitutive models only for olivine, but may be applied to other compositions (at your own risk). The main goal is to allow easy comparison between methods for calculating anelastic-dependent seismic properties, but the VBRc can also be used for calculating steady state viscosity, pure elastic (anharmonic) seismic properties and more. It can be used to fit and analyze experimental data, infer thermodynamic state from seismic measurements, predict measurable properties from geodynamic models, for example.
  </p>
</div>
```

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

``` [Phase mixing type]
Two phase :
* HS upper bound
* HS lower bound
* Modified Archie's law

Multiple phase :
* HS upper bound
* HS lower bound
* Generalized Archie's law
```
:::
