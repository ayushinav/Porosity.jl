abstract type AbstractAnelasticModel <: AbstractRockphyModel end

mutable struct RockphyAnelastic{T1, T2, T3, T4, T5, T6} <: AbstractRockphyResponse
    J1::T1
    J2::T2
    Qinv::T3
    M::T4
    V::T5
    Vave::T6
end

"""
    andrade_psp(T, P, dg, σ, ϕ, ρ, f)

Calculate anelastic properties stored in `RockPhyAnelastic` using the
Andrade model with pseudo-scaling per Jackson and Faul (2010)

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `dg`: Grain size (μm)
    - `σ` : Shear stress (GPa)
    - `ϕ` : Porosity
    - `ρ` : Density (kg/m³)
    - `f` : frequency

## Keyword Arguments

    - `params` : Various coefficients required for calculation.
    Also holds coefficients and the type of `RockphyElastic` model to be used.

    To investigate coefficients, call `default_params(Val{andrade_psp}())`. 
    To modify coefficients, check the relevant documentation page. This
    will also users to pick any particular type of `RockphyElastic` model, defaults to `anharmonic`.

## Usage

!!! note


**Make sure that the dimension of vector `f` is one more than the other parameters.
Check relevant tutorials. Note the transpose on `f` when making the model in the following eg.**

```julia
T = collect(1073.0f0:30:1373.0f0)
P = 2 .+ zero(T)
dg = collect(3.0f0:4.0f-1:7.0f0)
σ = collect(7.5f0:0.5f0:12.5f0) .* 1.0f-3
ϕ = collect(1.0f-2:1.0f-3:2.0f-2)
ρ = collect(3300.0f0:100.0f0:4300.0f0)

f = [1.0f0] #10f0 .^ collect(-10:1:0)

model = andrade_psp(T, P, dg, σ, ϕ, ρ, Ch2o_ol, f')

forward(model, [])
```

## References

  - Jackson and Faul, 2010, "Grainsize-sensitive viscoelastic relaxation in olivine: Towards a robust laboratory-based model for seismological application",
    Phys. Earth Planet. Inter.,
    https://doi.org/10.1016/j.pepi.2010.09.005
"""
mutable struct andrade_psp{T1, T2, T3, T4, T5, T6, T7} <: AbstractAnelasticModel
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    f::T7
end

"""
    eburgers_psp(T, P, dg, σ, ϕ, ρ, Ch2o_ol, T_solidus, f)

Calculate anelastic properties stored in `RockphyAnelastic` using the
Extended Burgers model with pseudo-scaling per Jackson and Faul (2010)

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `dg`: Grain size (μm)
    - `σ` : Shear stress (GPa)
    - `ϕ` : Porosity
    - `ρ` : Density (kg/m³)
    - `Ch2o_ol` : water concentration in olivine (in ppm)
    - `T_solidus` : Solidus temperature (K), only used when using `xfit_premelt` for viscosity calculations
    - `f` : frequency

## Keyword Arguments

    - `params` : Various coefficients required for calculation.
    Also holds coefficients and the type of `RockphyElastic` model and `RockphyViscous model` to be used.

    To investigate coefficients, call `default_params(Val{eburgers_psp}())`. 
    To modify coefficients, check the relevant documentation page. This
    will also users to pick any particular type of `RockphyElastic` model, defaults to `anharmonic`,
    as well as `RockphyViscous` model (for diffusion-derived viscosity), defaults to `xfit_premelt`

    `params` for `eburgers_psp` holds a few important fields:
        - `params_btype` : fitting parameters from Jackson and Faul (2010) to be used, defaults to `bg_only`. 
        Available options are : 
            - `bg_only` : multiple sample best high-temp background only fit
            - `bg_peak` : multiple sample best high-temp background + peak fit
            - `s6585_bg_only` : single sample 6585 fit, HTB only
            - `s6585_bg_peak` : single sample 6585 fit, HTB + dissipation peak

        - `melt_enhancement` : TODO

        -  `JF10_visc` : Whether to use scaling from Jackson and Faul (2010) for maxwell time calculations,
        otherwise calculate them using the `RockphyViscous` model provide. **Defaults to `true`.**

        - `integration_params` : Tells which integration option to be used, 
        and the number of points (Should not be touched ideally!)
        Available options are `quadgk`, `trapezoidal` and `simpson`, defaults to `quadgk`.

## Usage

!!! note


**Make sure that the dimension of vector `f` is one more than the other parameters.
Check relevant tutorials. Note the transpose on `f` when making the model in the following eg.**

```julia
T = collect(1073.0f0:30:1373.0f0)
P = 2 .+ zero(T)
dg = collect(3.0f0:4.0f-1:7.0f0)
σ = collect(7.5f0:0.5f0:12.5f0) .* 1.0f-3
ϕ = collect(1.0f-2:1.0f-3:2.0f-2)
T_solidus = 1473 .+ zero(T)
ρ = collect(3300.0f0:100.0f0:4300.0f0)
Ch2o_ol = zero(T)

f = [1.0f0] #10f0 .^ collect(-10:1:0)

model = eburgers_psp(T, P, dg, σ, ϕ, ρ, Ch2o_ol, T_solidus, f')

forward(model, [])
```

## References

  - Faul and Jackson, 2015, "Transient Creep and Strain Energy Dissipation: An Experimental Perspective",
    Ann. Rev. of Earth and Planetary Sci.,
    https://doi.org/10.1146/annurev-earth-060313-054732

  - Jackson and Faul, 2010, "Grainsize-sensitive viscoelastic relaxation in olivine:
    Towards a robust laboratory-based model for seismological application", Phys. Earth Planet. Inter.,
    https://doi.org/10.1016/j.pepi.2010.09.005
"""
mutable struct eburgers_psp{T1, T2, T3, T4, T5, T6, T7, T8, T9} <: AbstractAnelasticModel
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    Ch2o_ol::T7
    T_solidus::T8
    f::T9
end

eburgers_psp(T, P, dg, σ, ϕ, ρ, f) = eburgers_psp(T, P, dg, σ, ϕ, ρ, 0.0f0, 0.0f0, f) # TODO : args... ?

"""
    premelt_anelastic(T, P, dg, σ, ϕ, ρ, T_solidus, Ch2o_ol, f)

Calculate anelastic properties stored in `RockPhyAnelastic` using the
Master curve maxwell scaling per near-solidus parametrization of Yamauchi and Takei (2016),
with optional extension to include direct melt effects of Yamauchi and Takei (2024)

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `dg`: Grain size (μm)
    - `σ` : Shear stress (GPa)
    - `ϕ` : Porosity
    - `ρ` : Density (kg/m³)
    - `Ch2o_ol` : water concentration in olivine (in ppm)
    - `T_solidus` : Solidus temperature (K)
    - `f` : frequency

## Keyword Arguments

    - `params` : Various coefficients required for calculation.
    Also holds coefficients and the type of `RockphyElastic` model and `RockphyViscous model` to be used.

    To investigate coefficients, call `default_params(Val{xfit_premelt}())`. 
    To modify coefficients, check the relevant documentation page. This
    will also users to pick any particular type of `RockphyElastic` model, defaults to `anharmonic`.

    `params` for `premelt_anelastic` holds another important field:
        - `include_direct_melt_effect` : Whether to include the melt effect of Yamauchi and Takei (2024), defaults to false

## Usage

!!! note


**Make sure that the dimension of vector `f` is one more than the other parameters.
Check relevant tutorials. Note the transpose on `f` when making the model in the following eg.**

```julia
T = collect(1073.0f0:30:1373.0f0)
P = 2 .+ zero(T)
dg = collect(3.0f0:4.0f-1:7.0f0)
σ = collect(7.5f0:0.5f0:12.5f0) .* 1.0f-3
ϕ = collect(1.0f-2:1.0f-3:2.0f-2)
T_solidus = 1473 .+ zero(T)
ρ = collect(3300.0f0:100.0f0:4300.0f0)
Ch2o_ol = zero(T)

f = [1.0f0] #10f0 .^ collect(-10:1:0)

model = premelt_anelastic(T, P, dg, σ, ϕ, ρ, Ch2o_ol, T_solidus, f')

forward(model, [])
```

## References

  - Yamauchi and Takei, 2016, "Polycrystal anelasticity at near-solidus temperatures",
    J. Geophys. Res. Solid Earth,
    https://doi.org/10.1002/2016JB013316

  - Yamauchi and Takei, 2024, "Effect of Melt on Polycrystal Anelasticity",
    J. Geophys. Res. Solid Earth,
    https://doi.org/10.1029/2023JB027738
"""
mutable struct premelt_anelastic{T1, T2, T3, T4, T5, T6, T7, T8, T9} <:
               AbstractAnelasticModel
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    Ch2o_ol::T7
    T_solidus::T8
    f::T9
end

function premelt_anelastic(T, P, dg, σ, ϕ, ρ, T_solidus, f)
    premelt_anelastic(T, P, dg, σ, ϕ, ρ, 0.0f0, T_solidus, f)
end # TODO : args... ?

"""
    xfit_mxw(T, P, dg, σ, ϕ, ρ, T_solidus, Ch2o_ol, f)

Calculate anelastic properties stored in `RockPhyAnelastic` using the
Master curve maxwell scaling per McCarthy, Takei and Hiraga (2011)

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `dg`: Grain size (μm)
    - `σ` : Shear stress (GPa)
    - `ϕ` : Porosity
    - `ρ` : Density (kg/m³)
    - `Ch2o_ol` : water concentration in olivine (in ppm)
    - `T_solidus` : Solidus temperature (K), only used when using `xfit_premelt` for viscosity calculations
    - `f` : frequency

## Keyword Arguments

    - `params` : Various coefficients required for calculation. 
    Available options are `fit1` and `fit2`, defaults to `fit1`, i.e, `params_xfit_mxw.fit1`
    Also holds coefficients and the type of `RockphyElastic` model and `RockphyViscous model` to be used.

    To investigate coefficients, call `default_params(Val{xfit_premelt}())`. 
    To modify coefficients, check the relevant documentation page. This
    will also users to pick any particular type of `RockphyElastic` model, defaults to `anharmonic`,
    as well as `RockphyViscous` model (for diffusion-derived viscosity), defaults to `xfit_mxw`

## Usage

!!! note


**Make sure that the dimension of vector `f` is one more than the other parameters.
Check relevant tutorials. Note the transpose on `f` when making the model in the following eg.**

```julia
T = collect(1073.0f0:30:1373.0f0)
P = 2 .+ zero(T)
dg = collect(3.0f0:4.0f-1:7.0f0)
σ = collect(7.5f0:0.5f0:12.5f0) .* 1.0f-3
ϕ = collect(1.0f-2:1.0f-3:2.0f-2)
T_solidus = 1473 .+ zero(T)
ρ = collect(3300.0f0:100.0f0:4300.0f0)
Ch2o_ol = zero(T)

f = [1.0f0] #10f0 .^ collect(-10:1:0)

model = xfit_mxw(T, P, dg, σ, ϕ, ρ, Ch2o_ol, T_solidus, f')

forward(model, [])
```

## References

  - McCarthy, Takei, Hiraga, 2011, "Experimental study of attenuation and dispersion over a broad frequency range:

 2. The universal scaling of polycrystalline materials", Journal of Geophy Research,
    http://dx.doi.org/10.1029/2011JB008384
"""
mutable struct xfit_mxw{T1, T2, T3, T4, T5, T6, T7, T8, T9} <: AbstractAnelasticModel
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    Ch2o_ol::T7
    T_solidus::T8
    f::T9
end

"""
    andrade_analytical(T, P, dg, σ, ϕ, ρ, f)

Calculate anelastic properties stored in `RockPhyAnelastic` using the
Master curve maxwell scaling per McCarthy, Takei and Hiraga (2011)

## Arguments

    - `T` : Temperature of the rock (K)
    - `P` : Pressure (GPa)
    - `dg`: Grain size (μm)
    - `σ` : Shear stress (GPa)
    - `ϕ` : Porosity
    - `ρ` : Density (kg/m³)
    - `Ch2o_ol` : water concentration in olivine (in ppm), only used when using `HK2003` for viscosity calculations
    - `T_solidus` : Solidus temperature (K), only used when using `xfit_premelt` for viscosity calculations
    - `f` : frequency

## Keyword Arguments

    - `params` : Various coefficients required for calculation. 
    Also holds coefficients and the type of `RockphyElastic` model and `RockphyViscous model` to be used.

    To investigate coefficients, call `default_params(Val{xfit_premelt}())`. 
    To modify coefficients, check the relevant documentation page. This
    will also users to pick any particular type of `RockphyElastic` model, defaults to `anharmonic`,
    as well as `RockphyViscous` model (for diffusion-derived viscosity), defaults to `HK2003`

## Usage

!!! note


**Make sure that the dimension of vector `f` is one more than the other parameters.
Check relevant tutorials. Note the transpose on `f` when making the model in the following eg.**

```julia
T = collect(1073.0f0:30:1373.0f0)
P = 2 .+ zero(T)
dg = collect(3.0f0:4.0f-1:7.0f0)
σ = collect(7.5f0:0.5f0:12.5f0) .* 1.0f-3
ϕ = collect(1.0f-2:1.0f-3:2.0f-2)
ρ = collect(3300.0f0:100.0f0:4300.0f0)

f = [1.0f0] #10f0 .^ collect(-10:1:0)

model = andrade_analytical(T, P, dg, σ, ϕ, ρ, f')

forward(model, [])
```

## References

  - Andrade, 1910, "On the viscous flow in metals, and allied phenomena",
    Proceedings of the Royal Society of London,
    https://doi.org/10.1098/rspa.1910.0050

  - Cooper, 2002, "Seismic Wave Attenuation: Energy Dissipation in Viscoelastic Crystalline Solids",
    Reviews in mineralogy and geochemistry,
    https://doi.org/10.2138/gsrmg.51.1.253,
  - Lau and Holtzman, 2019, "“Measures of Dissipation in Viscoelastic Media” Extended:
    Toward Continuous Characterization Across Very Broad Geophysical Time Scales",
    Geophysical Research Letters,
    https://doi.org/10.1029/2019GL083529
"""
mutable struct andrade_analytical{T1, T2, T3, T4, T5, T6, T7, T8, T9} <:
               AbstractAnelasticModel
    T::T1
    P::T2
    dg::T3
    σ::T4
    ϕ::T5
    ρ::T6
    Ch2o_ol::T7
    T_solidus::T8
    f::T9
end

default_params(::Type{T}) where {T <: andrade_psp} = default_params_andrade_psp
default_params(::Type{T}) where {T <: eburgers_psp} = default_params_eburgers_psp
default_params(::Type{T}) where {T <: premelt_anelastic} = default_params_premelt_anelastic
default_params(::Type{T}) where {T <: xfit_mxw} = default_params_xfit_mxw
function default_params(::Type{T}) where {T <: andrade_analytical}
    default_params_andrade_analytical
end
