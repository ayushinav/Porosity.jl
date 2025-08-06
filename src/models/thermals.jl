"""
    solidus_Hirschmann2000(ps_nt)

returns pressure-dependent dry solidus (in K)

## Arguments

  - `P` : Pressure in GPa

## Optional Arguments

  - `params` : variables that parameterize the solidus

## Usage

The inputs should be in the form of a `NamedTuple`

```julia
P = [3.0, 4.0, 5.0]
ps_nt = (; P)
solidus_Hirschmann2000(ps_nt)
```

## References
"""
function solidus_Hirschmann2000(ps_nt)
    @unpack P = ps_nt
    T_solidus = @. 273 + 1108.08f0 + 139.44f0 * P - 5.904f0 * P * P
    return (; T_solidus)
end

"""
    solidus_Katz2003(ps_nt)

returns pressure-dependent dry solidus (in K)

## Arguments

  - `P` : Pressure in GPa

## Usage

The inputs should be in the form of a `NamedTuple`

```julia
P = [3.0, 4.0, 5.0]
ps_nt = (; P)
solidus_Katz2003(ps_nt)
```

## References

Katz, R. F., Spiegelman, M., & Langmuir, C. H. (2003) :
A new parameterization of hydrous mantle melting. Geochemistry, Geophysics, Geosystems, 4(9), 1–19.
https://doi.org/10.1029/2002GC000433
"""
function solidus_Katz2003(ps_nt)
    @unpack P = ps_nt
    T_solidus = @. 273 + 1085.7f0 + 132.9f0 * P - 5.1f0 * P * P
    return (; T_solidus)
end

"""
    ΔT_h2o_Katz2003(ps_nt)

updates solidus with the depressed one because of water in melt

## Arguments

  - `P` : Pressure (GPa)
  - `T_solidus` : solidus temperature (K)
  - `Ch2o_m` : melt fraction

# Usage

The inputs should be in the form of a `NamedTuple`

```julia
P = [3.0]
T_solidus = [1400.0] .+ 273
Ch2o_m = 0:100:10000
ps_nt = (; P, T_solidus, Ch2o_m)
ΔT_h2o_Katz2003(ps_nt)
```

## References

Katz, R. F., Spiegelman, M., & Langmuir, C. H. (2003) :
A new parameterization of hydrous mantle melting. Geochemistry, Geophysics, Geosystems, 4(9), 1–19.
https://doi.org/10.1029/2002GC000433
"""
function ΔT_h2o_Katz2003(ps_nt)
    @unpack Ch2o_m, P, T_solidus = ps_nt
    γ = 0.75f0
    K = 43.0f0

    Ch2o_m_sat = @. 12 * P^(0.6f0) + P
    # C_w = (Ch2o_m <= Ch2o_m_sat) ? (Ch2o_m) : Ch2o_m # worth investigating
    Ch2o_m_ = @. Ch2o_m * 1.0f-4
    dT = @. K * (Ch2o_m_)^γ
    T_solidus = @. T_solidus - dT
    return (; T_solidus)
end

function Dasgupta2007_core(Cco2_m_)
    if Cco2_m_ <= 25
        dT = 27.04f0 * Cco2_m_ + 1490.75f0 * log((100.0f0 - 1.18f0 * Cco2_m_) / 100.0f0)
    elseif Cco2_m_ > 25 && Cco2_m_ <= 37
        dTmax = 27.04f0 * 25.0f0 + 1490.75f0 * log((100.0f0 - 1.18f0 * 25.0f0) / 100.0f0)
        dT = dTmax + (160.0f0 - dTmax) / (37.0f0 - 25.0f0) * (Cco2_m_ - 25.0f0)
    elseif Cco2_m_ > 37
        dTmax = 27.04f0 * 25.0f0 + 1490.75f0 * log((100.0f0 - 1.18f0 * 25.0f0) / 100.0f0)
        dTmax = dTmax + (160.0f0 - dTmax)
        dT = dTmax + 150.0f0
    end
    return dT
end

"""
    ΔT_co2_Dasgupta2007(ps_nt)

updates solidus with the depressed one because of CO₂ in melt

## Arguments

  - `Cco2_m` : CO₂ conc. in melt (ppm.)
  - `T_solidus` : Temperature of solidus (K)

## Usage

The inputs should be in the form of a `NamedTuple`

```julia
P = [3.0]
T = [1400.0] .+ 273
Cco2_m = 0:100:10000
ps_nt = (; P, T, Cco2_m)
ΔT_co2_Dasgupta2007(ps_nt)
```

## References
"""
function ΔT_co2_Dasgupta2007(ps_nt)
    @unpack Cco2_m, T_solidus = ps_nt

    Cco2_m_ = @. Cco2_m * 1.0f-4 # wt %
    dT = @. Dasgupta2007_core(Cco2_m_)

    T_solidus = @. T_solidus - dT
    return (; T_solidus)
end

"""
    ΔT_co2_Dasgupta2013(ps_nt)

updates solidus with the depressed one because of CO₂ in melt

## Arguments

  - `Cco2_m` : CO₂ conc. in melt (ppm.)
  - `T_solidus` : Temperature of solidus (K)

## Usage

The inputs should be in the form of a `NamedTuple`

```julia
P = [3.0]
T = [1400.0] .+ 273
Cco2_m = 0:100:10000
ps_nt = (; P, T, Cco2_m)
ΔT_co2_Blatter2022(ps_nt)
```

## References

Dasgupta, R., Mallik, A., Tsuno, K. et al. (2013) :
Carbon-dioxide-rich silicate melt in the Earth’s upper mantle.
Nature 493, 211–215 (2013). https://doi.org/10.1038/nature11731
"""
function ΔT_co2_Dasgupta2013(ps_nt)
    @unpack Cco2_m, T_solidus = ps_nt

    Cco2_m_ = @.Cco2_m * 1.0f-4

    #=
    2 GPa: a = 19.21; b = 1491.37; c = 0.86 (this study) 
    3 GPa: a = 27.04; b = 1490.75; c = 1.18 (ref. 5) 
    4 GPa: a = 31.90; b = 1469.92; c = 1.31 (this study) 
    5 GPa: a = -5.01; b = 1514.84; c = -1.23 (this study) 
    =#

    a = 19.21f0
    b = 1491.37f0
    c = 0.86f0

    ΔT_co2 = @. a * Cco2_m_ + b * log(1 - c * Cco2_m_ * 1.0f-2)
    T_co2 = @. T_solidus - ΔT_co2
    return (; T_solidus=T_co2)
end

"""
    ΔT_h2o_Blatter2022(ps_nt)

updates solidus with the depressed one because of water in melt

## Arguments

  - `Ch2o_m` : water conc. in melt (ppm.)
  - `T_solidus` : Temperature of solidus (K)

## Usage

The inputs should be in the form of a `NamedTuple`

```julia
P = [3.0]
T = [1400.0] .+ 273
Ch2o_m = 0:100:10000
ps_nt = (; P, T, Ch2o_m)
ΔT_h2o_Blatter2022(ps_nt)
```

## References

Blatter D, Naif S, Key K, Ray A (2022) :
A plume origin for hydrous melt at the lithosphere-asthenosphere boundary.
Nature. 2022 Apr;604(7906):491-494. doi: 10.1038/s41586-022-04483-w
"""
function ΔT_h2o_Blatter2022(ps_nt)
    @unpack Ch2o_m, T_solidus = ps_nt

    # D = 0.005
    Ch2o_m_ = @. Ch2o_m * 1.0f-4

    M = 59.0f0
    ΔS = 0.4

    X_OH = @. 2M * Ch2o_m_ / 100 / (18.02 + Ch2o_m_ / 100 * (2M - 18.02))
    T_wet = @. 273 +
               (T_solidus - 273) * inv(1 - MT.gas_R * 1.0f3 / (M * ΔS) * log(1 - X_OH))

    return (; T_solidus=T_wet)
end

function get_Cco2_m_core(ϕ, Cco2, Cco2_sat)
    Cco2_m = Cco2 * inv(ϕ)
    if Cco2_m > Cco2_sat
        return Cco2_sat
    else
        return Cco2_m
    end
end

"""
    get Cco2_m(ps_nt)

returns the amount of CO₂ in melt

## Arguments

  - `ϕ` : melt fraction
  - `Cco2` : bulk CO₂ conc.
  - `Cco2_sat` : saturation limit of CO₂ conc. in melt

## Usage

The inputs should be in the form of a `NamedTuple`

```julia
ϕ = 0.05
Cco2 = 1000
ps_nt = (; ϕ, Cco2)
get_Cco2_m(ps_nt)
```
"""
function get_Cco2_m(ps_nt)
    @unpack ϕ, Cco2, Cco2_sat = ps_nt
    Cco2_m = @. get_Cco2_m_core(ϕ, Cco2, Cco2_sat)
    return (; Cco2_m)
end

"""
    get Ch2o_m(ps_nt)

returns the amount of CO₂ in melt and water in the form of `NamedTuple`

## Arguments

  - `ϕ` : melt fraction
  - `Cco2` : bulk water conc.

## Usage

The inputs should be in the form of a `NamedTuple`

```julia
ϕ = 0.05
Ch2o = 1000
ps_nt = (; ϕ, Ch2o)
get_Cco2_m(ps_nt)
```
"""
function get_Ch2o_m(ps_nt)
    @unpack ϕ, Ch2o, D = ps_nt
    Ch2o_m = @. Ch2o * inv(D + ϕ * (1 - D))
    Ch2o_ol = @. 0 * Ch2o_m
    Ch2o_ol = @. D * Ch2o_m
    return (; Ch2o_m, Ch2o_ol)
end

function get_melt_fraction_core(
        Ch2o, Cco2, T_solidus, P, D, H2O_suppress_fn, CO2_suppress_fn)
    function f(u, p)
        Ch2o_m = get_Ch2o_m((; ϕ=u, p.Ch2o, p.D)).Ch2o_m
        Cco2_m = get_Cco2_m((; ϕ=u, p.Cco2)).Cco2_m

        T_new_H2O = H2O_suppress_fn((; p..., Ch2o_m)).T_solidus
        T_new_CO2 = H2O_suppress_fn((; p..., Cco2_m)).T_solidus
        ΔT = 2p.T_solidus - T_new_H2O - T_new_CO2
        dTdF = -40 * p.P + 450

        return u * dTdF / ΔT - 1
    end

    prob_init = IntervalNonlinearProblem(
        f, (1.0f-15, 1.0f0), (; Ch2o, Cco2, T_solidus, P, D))
    sol = solve(prob_init)
    return sol.u
end

"""
    get_melt_fraction(ps_nt; H2O_suppress_fn = ΔT_h2o_Blatter2022, CO2_suppress_fn= ΔT_co2_Blatter2022)

returns the melt fraction that is thermodynamically stable at given water conc. and CO₂ conc. for a given solidus temperature and pressure

## Arguments

  - `Cco2` : bulk CO₂ conc. (ppm)
  - `Ch2o` : bulk water conc. (ppm)
  - `T_solidus` : solidus temperature (K)
  - `P` : Pressure (GPa)
  - `D` : water partition coefficient

## Keyword Arguments

  - `H2O_suppress_fn` : function to calculate suppressed solidus due to presence of water, defaults to `ΔT_h2o_Blatter2022`
  - `CO2_suppress_fn` : function to calculate suppressed solidus due to presence of CO₂, defaults to `ΔT_co2_Dasgupta2013`

## Usage

The inputs should be in the form of a `NamedTuple`

```julia
P = [3.0]
T_solidus = [1400.0] .+ 273
Ch2o = 0:100:10000
Cco2 = 1000
ps_nt = (; P, T_solidus, Ch2o_m, Cco2_m)
get_melt_fraction(ps_nt)
```
"""
function get_melt_fraction(
        ps_nt; H2O_suppress_fn=ΔT_h2o_Blatter2022, CO2_suppress_fn=ΔT_co2_Dasgupta2013)
    @unpack Cco2, Ch2o, T_solidus, P, D = ps_nt

    if Ch2o ∈ keys(ps_nt)
        Ch2o = ps_nt.Ch2o
    else
        Ch2o = 0.0f0
    end

    if Cco2 ∈ keys(ps_nt)
        Cco2 = ps_nt.Cco2
    else
        Cco2 = 0.0f0
    end

    ϕ = broadcast(get_melt_fraction_core, Ch2o, Cco2, T_solidus,
        P, D, H2O_suppress_fn, CO2_suppress_fn)

    return (; ϕ)
end
