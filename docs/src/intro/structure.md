# Structure

One of the goals of this package is to provide an easy and intiutive way to compute the geophysical observables. In order to that, it is helpful to understand the structure of codes. The following serves as a cheatsheet to understand the structure.

## Building blocks

We begin with the example of [SEO3](@ref), a simple model that estimates the temperature dependent conductivity of olivine. It is a `AbstractCondModel`, which in simple terms groups all the models that estimate conductivity. Similarly, we have `AbstractElasticModel`, `AbstractViscousModel` and `AbstractAnelasticModel`. These are the main building blocks and are shown at the top in the following flowchart.

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

## Mixing phases
More often than not, we want to compute the bulk conductivity in the presence of fluids. This can be accompolished with [`two_phase_model`](@ref) and [`multi_phase_model`](@ref).
[Relevant tutorial](../tutorials/mixing_phases.md)

## Multi-rock physics
A further goal of the package is to encourage joint stochastic inference of more than one data types. For that purpose, we provide [`multi_rp_model`](@ref).
[Relevant tutorial](../tutorials/combine_models.md)

## Tuning rock physics
It happens in a lot of cases that one of the inputs to a rock physics model is dependent on the other, eg. melt fraction depends on temperature of the rock and its solidus temperature. In those cases where we want to tune one of the models with respect to others, [`tune_rp_model`] will come handy.
[Relevant tutorial](../tutorials/tune_rp.md)
