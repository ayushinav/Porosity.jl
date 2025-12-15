"""
    solidus_Hirschmann2000(ps_nt)

returns pressure-dependent dry solidus (in K)

## Arguments

  - `P` : Pressure in GPa

## Optional Arguments

  - `params` : variables that parameterize the solidus

## Usage

The inputs should be in the form of a `NamedTuple`, and the output will be a `NamedTuple` as well with `T_solidus` as the only field

```jldoctest; output = false
P = [3.0, 4.0, 5.0]
ps_nt = (; P)
solidus_Hirschmann2000(ps_nt)

# output

(T_solidus = [1746.2639651298523, 1844.3759689331055, 1930.6799731254578],)
```

## References

Marc M. Hirschmann (2000) : Mantle solidus: Experimental constraints and the effects of peridotite composition
https://doi.org/10.1029/2000GC000070
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

The inputs should be in the form of a `NamedTuple`, and the output will be a `NamedTuple` as well with `T_solidus` as the only field

```jldoctest; output = false
P = [3.0, 4.0, 5.0]
ps_nt = (; P)
solidus_Katz2003(ps_nt)

# output

(T_solidus = [1711.499933719635, 1808.6999282836914, 1895.6999230384827],)
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

## Usage

The inputs should be in the form of a `NamedTuple`, and the output will be a `NamedTuple` as well with `T_solidus` as the only field

```jldoctest; output = false
P = [3.0]
T_solidus = [1400.0] .+ 273
Ch2o_m = 0:200:1000
ps_nt = (; P, T_solidus, Ch2o_m)
ΔT_h2o_Katz2003(ps_nt)

# output

(T_solidus = [1673.0, 1670.7131328582764, 1669.1539633274078, 1667.7870697975159, 1666.531762599945, 1665.3533987998962],)
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

    etype_ = eltype(P .+ Ch2o_m)

    Ch2o_m_sat = @. etype_(12 * P^(0.6f0) + P)
    Ch2o_m_ = @. ifelse(Ch2o_m * 1.0f-4 < Ch2o_m_sat, etype_(Ch2o_m * 1.0f-4), etype_(Ch2o_m_sat))
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

```jldoctest; output = false
P = [3.0]
T_solidus = [1400.0] .+ 273
Cco2_m = 0:200:1000
ps_nt = (; P, T_solidus, Cco2_m)
ΔT_co2_Dasgupta2007(ps_nt)

# output

(T_solidus = [1673.0, 1672.81102091074, 1672.6222137212753, 1672.4334008693695, 1672.2446709871292, 1672.0560241937637],)
```

## References

Rajdeep Dasgupta; Marc M. Hirschmann; Neil D. Smith (2007) :
Water follows carbon: CO 2 incites deep silicate melting and dehydration beneath mid-ocean ridges, Geology, 2007
https://doi.org/10.1130/G22856A.1
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

```jldoctest; output = false
P = [3.0]
T_solidus = [1400.0] .+ 273
Cco2_m = 0:200:1000
ps_nt = (; P, T_solidus, Cco2_m)
ΔT_co2_Dasgupta2013(ps_nt)

# output

(T_solidus = [1673.0, 1672.8723660707474, 1672.7446873784065, 1672.6171417236328, 1672.4896402359009, 1672.3620940446854],)
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

```jldoctest; output = false
P = [3.0]
T_solidus = [1400.0] .+ 273
Ch2o_m = 0:200:1000
ps_nt = (; P, T_solidus, Ch2o_m)
ΔT_h2o_Blatter2022(ps_nt)

# output

(T_solidus = [1673.0, 1672.2288240004204, 1671.459057866923, 1670.6906970192301, 1669.9237371410397, 1669.1581736299327],)
```

## References

Blatter D, Naif S, Key K, Ray A (2022) :
A plume origin for hydrous melt at the lithosphere-asthenosphere boundary.
Nature. 2022 Apr;604(7906):491-494. doi: 10.1038/s41586-022-04483-w
"""
function ΔT_h2o_Blatter2022(ps_nt)
    @unpack Ch2o_m, T_solidus = ps_nt

    # D = 0.005
    Ch2o_m_ = @. Ch2o_m * 1.0f-6

    M = 59.0f0
    ΔS = 0.4

    X_OH = @. 2M * Ch2o_m_ / (18.02 + Ch2o_m_ * (2M - 18.02))
    T_wet = @. (T_solidus) * inv(1 - gas_R * 1.0f3 / (M * ΔS) * log(1 - X_OH))

    return (; T_solidus=T_wet)
end

function get_Cco2_m_core(ϕ, Cco2, Cco2_sat)
    type_ = typeof(ϕ + Cco2 + Cco2_sat)
    if ϕ == 0
        return type_(0.0f0)
    elseif Cco2 * inv(ϕ) > Cco2_sat
        return type_(Cco2_sat)
    else
        return type_(Cco2 * inv(ϕ))
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

```jldoctest
ϕ = 0.05
Cco2 = 1000
Cco2_sat = 38e4
ps_nt = (; ϕ, Cco2, Cco2_sat)
get_Cco2_m(ps_nt)

# output

(Cco2_m = 20000.0,)
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
  - `D` : partition coefficient

## Usage

The inputs should be in the form of a `NamedTuple`

```jldoctest
ϕ = 0.05
Ch2o = 1000
D = 0.005
ps_nt = (; ϕ, Ch2o, D)
get_Ch2o_m(ps_nt)

# output

(Ch2o_m = 18264.840182648404, Ch2o_ol = 91.32420091324202)
```
"""
function get_Ch2o_m(ps_nt)
    @unpack ϕ, Ch2o, D = ps_nt
    Ch2o_m = @. Ch2o * inv(D + ϕ * (1 - D))
    Ch2o_ol = @. D * Ch2o_m
    return (; Ch2o_m, Ch2o_ol)
end

function f_melt(u, p, H2O_suppress_fn, CO2_suppress_fn)
    T, T_solidus, Ch2o, Cco2, Cco2_sat, P, D = p
    # ps = (; T, T_solidus, Ch2o, Cco2, Cco2_sat, P, D)

    Ch2o_m = Ch2o * inv(D + u * (1 - D))
    Cco2_m = get_Cco2_m((; ϕ=u, Cco2, Cco2_sat)).Cco2_m

    T_new_H2O = H2O_suppress_fn((; T, T_solidus, Ch2o, Cco2, Cco2_sat, P, D, Ch2o_m)).T_solidus
    T_new_CO2 = CO2_suppress_fn((; T, T_solidus, Ch2o, Cco2, Cco2_sat, P, D, Cco2_m)).T_solidus
    T_solidus_new = T_solidus - (2T_solidus - T_new_H2O - T_new_CO2)
    ΔT = max(zero(T - T_solidus_new), T - T_solidus_new)
    dTdF = -40 * P + 450
    return u - ΔT / dTdF
end

function get_melt_fraction_core(
        T, T_solidus, Ch2o, Cco2, Cco2_sat, P, D, H2O_suppress_fn, CO2_suppress_fn)
    f(u, p) = f_melt(u, p, H2O_suppress_fn, CO2_suppress_fn)

    p = [T, T_solidus, Ch2o, Cco2, Cco2_sat, P, D]

    f1 = f(1.0f-15, p)
    f2 = f(1.0f0, p)
    etype_ = eltype(p)

    if f1 * f2 > 0
        return zero(etype_)
    else
        prob_init = IntervalNonlinearProblem{false}(f, (
            exp10(-15 * one(etype_)), one(etype_)), p)
        sol = solve(prob_init)
        return etype_(sol.u)
    end
end

"""
    get_melt_fraction(ps_nt; H2O_suppress_fn = ΔT_h2o_Blatter2022, CO2_suppress_fn= ΔT_co2_Blatter2022)

returns the melt fraction that is thermodynamically stable at given water conc. and CO₂ conc. for a given dry solidus temperature and pressure

## Arguments

  - `Cco2` : bulk CO₂ conc. (ppm), defaults to 0 ppm
  - `Cco2_sat` : saturation limit of melt CO₂ conc. (ppm), defaults to 38 wt %
  - `Ch2o` : bulk water conc. (ppm), defaults to 0 ppm
  - `T` : temperature (K)
  - `T_solidus` : solidus temperature (K)
  - `P` : Pressure (GPa)
  - `D` : water partition coefficient between water and olivine, defaults to 0.005

## Keyword Arguments

  - `H2O_suppress_fn` : function to calculate suppressed solidus due to presence of water, defaults to `ΔT_h2o_Blatter2022`
  - `CO2_suppress_fn` : function to calculate suppressed solidus due to presence of CO₂, defaults to `ΔT_co2_Dasgupta2013`

## Usage

The inputs should be in the form of a `NamedTuple`

```jldoctest
P = [3.0]
T = (1440:10:1480) .+ 273.0
T_solidus = [1450.0] .+ 273
Ch2o = 100.0
Cco2 = 1000.0

ps_nt = (; P, T, T_solidus, Ch2o, Cco2)
get_melt_fraction(ps_nt)

# output

(ϕ = [0.041345720244193085, 0.05481061353960702, 0.072300827427471, 0.09337629457529527, 0.1171677298066955],)
```
"""
function get_melt_fraction(ps_nt; H2O_suppress_fn=ΔT_h2o_Blatter2022, CO2_suppress_fn=ΔT_co2_Dasgupta2013)
    ps_nt = (; Ch2o=0.0f0, Cco2=0.0f0, Cco2_sat=38.0f4, D=0.005f0, ps_nt...)
    @unpack Cco2, Cco2_sat, Ch2o, T, T_solidus, P, D = ps_nt

    ϕ = broadcast(get_melt_fraction_core, T, T_solidus, Ch2o, Cco2,
        Cco2_sat, P, D, H2O_suppress_fn, CO2_suppress_fn)

    return (; ϕ)
end
