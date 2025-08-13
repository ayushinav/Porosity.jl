const global boltz_k = 8.617f-5
const global charge_e = 1.602f-19
const global gas_R = 0.008314f0

# abstract type AbstractMineralModel end
# abstract type AbstractMeltModel end

abstract type AbstractCondModel <: AbstractRockphyModel end

# ========================================================================================================== 
# olivine

mutable struct RockphyCond{T} <: AbstractRockphyResponse
    σ::T
end

"""
    SEO3(T)

Electrical conductivity model for olivine dependent on temperature.

## Arguments

    - `T` : Temperature of olivine (in K)

## Usage

```julia
model = SEO3(1000 + 273.0)

log_cond = forward(model, [])
```

## References

  - Constable, S (2006), "SEO3: A new model of olivine electrical conductivity", Geophysical Journal International,
    Volume 166, Issue 1, July 2006, Pages 435–437, https://doi.org/10.1111/j.1365-246X.2006.03041.x
"""
mutable struct SEO3{F} <: AbstractCondModel
    T::F
end

"""
    UHO2014(T, Ch2o_ol)

Electrical conductivity model for olivine dependent on temperature and water concentration.

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_ol` : water concentration in olivine (in ppm)

## References

  - Gardés, E., F. Gaillard, and P. Tarits (2014), "Toward a unified hydrous olivine electrical conductivity law",
    Geochem. Geophys. Geosyst., 15, 4984–5000, doi:10.1002/2014GC005496.

## Usage

```julia
model = UHO2014(1000 + 273.0, 2e4)

log_cond = forward(model, [])
```
"""
mutable struct UHO2014{F1, F2} <: AbstractCondModel
    T::F1
    Ch2o_ol::F2
end

"""
    Jones2012(T, Ch2o_ol)

Electrical conductivity model for olivine dependent on temperature and water concentration.

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_ol` : water concentration in olivine (in ppm)

## Usage

```julia
model = Jones2012(1000 + 273.0, 2e4)

log_cond = forward(model, [])
```

## References

  - Jones, A. G., J. Fullea, R. L. Evans, and M. R. Muller (2012), "Water in cratonic lithosphere:
    Calibrating laboratory-determined models of electrical conductivity of mantle minerals using geophysical and petrological observations",
    Geochem. Geophys. Geosyst., 13, Q06010, doi:10.1029/2012GC004055.
"""
mutable struct Jones2012{F1, F2} <: AbstractCondModel
    T::F1
    Ch2o_ol::F2
end

"""
    Poe2010(T, Ch2o_ol)

Electrical conductivity model for olivine dependent on temperature and water concentration.

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_ol` : water concentration in olivine (in ppm)

## Usage

```julia
model = Poe2010(1000 + 273.0, 2e4)

log_cond = forward(model, [])
```

## References

  - Brent T. Poe, Claudia Romano, Fabrizio Nestola, Joseph R. Smyth (2010),
    "Electrical conductivity anisotropy of dry and hydrous olivine at 8GPa",
    Physics of the Earth and Planetary Interiors,Volume 181, Issues 3–4, 2010, Pages 103-111, ISSN 0031-9201,
    https://doi.org/10.1016/j.pepi.2010.05.003.
"""
mutable struct Poe2010{F1, F2} <: AbstractCondModel
    T::F1
    Ch2o_ol::F2
end

"""
    Wang2006(T, Ch2o_ol)

Electrical conductivity model for olivine dependent on temperature and water concentration.

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_ol` : water concentration in olivine (in ppm)

## Usage

model = Wang2006(1000 + 273., 2e4)

log_cond = forward(model, [])

## References

  - Wang, D., Mookherjee, M., Xu, Y. et al. (2006),
    "The effect of water on the electrical conductivity of olivine", Nature 443, 977–980 (2006),
    doi: https://doi.org/10.1038/nature05256
"""
mutable struct Wang2006{F1, F2} <: AbstractCondModel
    T::F1
    Ch2o_ol::F2
end

"""
    Yoshino2009(T, Ch2o_ol)

Electrical conductivity model for olivine dependent on temperature and water concentration.

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_ol` : water concentration in olivine (in ppm)

## Usage

```julia
model = Yoshino2009(1000 + 273.0, 2e4)

log_cond = forward(model, [])
```

## References

  - Takashi Yoshino, Takuya Matsuzaki, Anton Shatskiy, Tomoo Katsura (2009),
    "The effect of water on the electrical conductivity of olivine aggregates and its implications for the electrical structure of the upper mantle,
    Earth and Planetary Science Letters",
    Volume 288, Issues 1–2, 2009, Pages 291-300, ISSN 0012-821X,
    https://doi.org/10.1016/j.epsl.2009.09.032
"""
mutable struct Yoshino2009{F1, F2} <: AbstractCondModel
    T::F1
    Ch2o_ol::F2
end

"""
    const_matrix(σ)

Fixed electrical conductivity model for the phase.

## Arguments

  - `σ` : Conductivity of the phase

## Usage

```julia
model = const_matrix(1000.0)

log_cond = forward(model, [])
```
"""
mutable struct const_matrix{F} <: AbstractCondModel
    σ::F
end

# ========================================================================================================== 
# melts
"""
    Ni2011(T, Ch2o_m)

Electrical conductivity model for basaltic melt dependent on Temperature and water content in melt.

## Arguments

  - `T` : Temperature of melt (should be greater than 1146.8 K)
  - `Ch2o_m` : water concentration in melt (in ppm)

## Usage

```julia
model = Ni2011(1000 + 273.0, 2e4)

log_cond = forward(model, [])
```

## References

  - Ni, H., Keppler, H. & Behrens, H. (2011),
    "Electrical conductivity of hydrous basaltic melts: implications for partial melting in the upper mantle.",
    Contrib Mineral Petrol 162, 637–650 (2011), doi: https://doi.org/10.1007/s00410-011-0617-4
"""
mutable struct Ni2011{F1, F2} <: AbstractCondModel
    T::F1
    Ch2o_m::F2
    function Ni2011(T, Ch2o_m)
        if any(T .< params_Ni2011.T_corr)
            @warn "T (= $T K) should be greater than $(params_Ni2011.T_corr) K for `Ni2011` otherwise erroneous values are obtained"
        end
        return new{typeof(T), typeof(Ch2o_m)}(T, Ch2o_m)
    end
end

"""
    Sifre2014(T, Ch2o_m, Cco2_m)

Electrical conductivity model for melt dependent on Temperature, water content and CO₂ content in melt.

## Arguments

  - `T` : Temperature of melt (should be greater than 1146.8 K)
  - `Ch2o_m` : water concentration in melt (in ppm)
  - `Cco2_m` : Co2 concentration in melt (in ppm)

## Usage

```julia
model = Sifre2014(1000 + 273.0, 2e4, 2e4)

log_cond = forward(model, [])
```

## References

  - Sifré, D., Gardés, E., Massuyeau, M. et al. (2014), "Electrical conductivity during incipient melting in the oceanic low-velocity zone."
    Nature 509, 81–85 (2014), doi: https://doi.org/10.1038/nature13245
"""
mutable struct Sifre2014{F1, F2, F3} <: AbstractCondModel
    T::F1
    Ch2o_m::F2
    Cco2_m::F3
end

"""
    Gaillard2008(T)

Electrical conductivity model for melt dependent on Temperature.

## Usage

```julia
model = Gaillard2008(1000 + 273.0)

log_cond = forward(model, [])
```

## Arguments

  - `T` : Temperature of melt (should be greater than 1146.8 K)

## References

  - Gaillard, Fabrice & Malki, Mohammed & Iacono-Marziano, Giada & Pichavant, Michel & Scaillet, Bruno. (2008),
    "Carbonatite Melts and Electrical Conductivity in the Asthenosphere",
    Science (New York, N.Y.). 322. 1363-5, doi: 10.1126/science.1164446.
"""
mutable struct Gaillard2008{F1} <: AbstractCondModel
    T::F1
end

# ========================================================================================================== 
# orthopyroxene

"""
    Dai_Karato2009(T, Ch2o_opx)

Electrical conductivity model for olivine dependent on temperature and water concentration.

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_opx` : water concentration in orthopyroxene (in ppm)

## References

  - todo

## Usage

```julia
model = Dai_Karato2009(1000 + 273.0, 2e4)

log_cond = forward(model, [])
```
"""
mutable struct Dai_Karato2009{F1, F2} <: AbstractCondModel
    T::F1
    Ch2o_opx::F2
end

"""
    Zhang2012(T, Ch2o_opx)

Electrical conductivity model for olivine dependent on temperature and water concentration.

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_opx` : water concentration in olivine (in ppm)

## References

  - todo

## Usage

```julia
model = Zhang2012(1000 + 273.0, 2e4)

log_cond = forward(model, [])
```
"""
mutable struct Zhang2012{F1, F2} <: AbstractCondModel
    T::F1
    Ch2o_opx::F2
end

# ========================================================================================================== 
# clinopyroxene

"""
    Yang2011(T, Ch2o_cpx)

Electrical conductivity model for olivine dependent on temperature and water concentration.

## Arguments

  - `T` : Temperature of olivine (in K)
  - `Ch2o_cpx` : water concentration in olivine (in ppm)

## References

  - todo

## Usage

```julia
model = Yang2011(1000 + 273.0, 2e4)

log_cond = forward(model, [])
```
"""
mutable struct Yang2011{F1, F2} <: AbstractCondModel
    T::F1
    Ch2o_cpx::F2
end

default_params(::Type{T}) where {T <: SEO3} = default_params_SEO3
default_params(::Type{T}) where {T <: UHO2014} = default_params_UHO2014
default_params(::Type{T}) where {T <: Jones2012} = default_params_Jones2012
default_params(::Type{T}) where {T <: Poe2010} = default_params_Poe2010
default_params(::Type{T}) where {T <: Wang2006} = default_params_Wang2006
default_params(::Type{T}) where {T <: Yoshino2009} = default_params_Yoshino2009
default_params(::Type{T}) where {T <: Ni2011} = default_params_Ni2011
default_params(::Type{T}) where {T <: Sifre2014} = default_params_Sifre2014
default_params(::Type{T}) where {T <: Gaillard2008} = default_params_Gaillard2008
default_params(::Type{T}) where {T <: Dai_Karato2009} = default_params_Dai_Karato2009
default_params(::Type{T}) where {T <: Zhang2012} = default_params_Zhang2012
default_params(::Type{T}) where {T <: Yang2011} = default_params_Yang2011
default_params(::Type{T}) where {T <: const_matrix} = default_params_const_matrix
