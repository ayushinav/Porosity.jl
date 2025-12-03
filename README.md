# Porosity.jl

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://ayushinav.github.io/Porosity.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://ayushinav.github.io/Porosity.jl/dev/)

[![Coverage](https://codecov.io/gh/ayushinav/MT.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/ayushinav/MT.jl)
[![Coverage](https://coveralls.io/repos/github/ayushinav/MT.jl/badge.svg?branch=main)](https://coveralls.io/github/ayushinav/MT.jl?branch=main)

`Porosity.jl` is a rock physics modeling and inference library written in julia, providing performant, AD compatible and scalable codes.

Available methods are : 

## Conductivity models

### olivine

* Jones2012
* Poe2010
* SEO3
* UHO2014
* Wang2006
* Yang2011
* Yoshino2009

### melt
* Gaillard2008
* Ni2011
* Sifre2014

### orthopyroxene

* Dai_Karato2009
* Zhang2012

## Elastic models

* anharmonic
* anhormonic_poro
* SLB2005

## Viscosity models

* HK2003
* HZK2011
* xfit_premelt

```mermaid
graph TD;
    subgraph L[ ]
        AbstractCondModel
        AbstractElasticModel
        AbstractViscousModel
        AbstractAnelasticModel
    end

    AbstractCondModel --> Matrix_conductivity;
    AbstractCondModel --> Melt_conductivity;
    AbstractCondModel -.-> m1[" "];
    AbstractCondModel -.-> m2[" "];
    AbstractCondModel -.-> multi_phase_model;
    m1 -.-> multi_phase_model
    m2 -.-> multi_phase_model
    multi_phase_model --> combine_rp_model
    Matrix_conductivity --> two_phase_model;
    Melt_conductivity --> two_phase_model;
    two_phase_model --> combine_rp_model;
    AbstractCondModel --> combine_rp_model;
    AbstractElasticModel --> combine_rp_model;
    AbstractViscousModel --> combine_rp_model;
    AbstractAnelasticModel --> combine_rp_model;
    combine_rp_model --> tune_rp_model;

    style L fill:#fff,stroke:#fff
```
