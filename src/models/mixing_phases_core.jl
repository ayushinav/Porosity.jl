abstract type phase_mixing end

"""
    HS1962_plus

Hashim-Strikman upper bound for mixing 2 phases

## Usage

```julia
m = two_phase_modelType(Yoshino2009, Sifre2014, HS1962_plus())
```

Check out the relevant rock physics documentation.

## References

  - Paul W. J. Glover (2010), "A generalized Archie's law for n phases",
    Geophysics 2010; 75 (6): E247–E265. doi: https://doi.org/10.1190/1.3509781
"""
mutable struct HS1962_plus <: phase_mixing end

"""
    HS1962_minus

Hashim-Strikman lower bound for mixing 2 phases

## Usage

```julia
m = two_phase_modelType(Yoshino2009, Sifre2014, HS1962_minus())
```

Check out the relevant rock physics documentation.

## References

  - Paul W. J. Glover (2010), "A generalized Archie's law for n phases",
    Geophysics 2010; 75 (6): E247–E265, doi: https://doi.org/10.1190/1.3509781
"""
mutable struct HS1962_minus <: phase_mixing end

"""
    MAL

Modified Archie's law for mixing 2 phases

## Usage

```julia
m = two_phase_modelType(Yoshino2009, Sifre2014, MAL(0.2))
```

Check out the relevant rock physics documentation.

## References

    - Glover, P. W. J., Hole, M. J., & Pous, J. (2000),
    "A Modified Archie’s Law for two conducting phases", Earth and Planetary Science Letters, 180(3–4), 369–383, 
    doi: https://doi.org/10.1016/S0012-821X(00)00168-0
"""
mutable struct MAL{T} <: phase_mixing
    m_MAL::T
end

"""
Hashim-Strikman upper bound for multiple phases

## References

    - Berryman, J. G. (1995), Mixture Theories for Rock Properties
     https://doi.org/10.1029/RF003p0205
"""
mutable struct HS_plus_multi_phase <: phase_mixing end

"""
Hashim-Strikman lower bound for multiple phases

## References

    - Berryman, J. G. (1995), Mixture Theories for Rock Properties
     https://doi.org/10.1029/RF003p0205
"""
mutable struct HS_minus_multi_phase <: phase_mixing end

"""
    GAL(m)

Generalized Archie's law for multiple phases

## Usage

## References

    - Glover, P.W.J (2010), A generalized Archie’s law for n phases
    https://doi.org/10.1190/1.3509781
"""
mutable struct GAL{T} <: phase_mixing
    m_GAL::T
end

"""
    single_phase

Single phase only conductivity. Assumes the rock matrix is composed of a single phase only.
"""
mutable struct single_phase <: phase_mixing end

# mixing functions

two_phase_mix_types = Union{HS1962_plus, HS1962_minus, MAL}
multi_phase_mix_types = Union{HS_plus_multi_phase, HS_minus_multi_phase, GAL}

function mix_models(σs, ϕ, ::HS1962_plus)
    σ_max = maximum(σs)
    σ_min = minimum(σs)
    phi = first(ϕ)

    num = 3 * (1 - phi) * (σ_max - σ_min) # numerator
    den = 3 * σ_max - phi * (σ_max - σ_min) # denominator
    esig = σ_max * (1 - (num / den))

    return esig
end

function mix_models(σs, ϕ, ::HS1962_minus)
    σ_max = maximum(σs)
    σ_min = minimum(σs)
    phi = first(ϕ)

    num = 3 * (phi) * (σ_max - σ_min) # numerator
    den = 3 * σ_min + (1 - phi) * (σ_max - σ_min) # denominator
    esig = σ_min * (1 + (num / den))

    return esig
end

function mix_models(σs, ϕ, mal::MAL)
    σ_fluid = (σs[2])
    σ_matrix = (σs[1])

    phi = first(ϕ)
    sig = σ_fluid

    if phi < 1
        p = log10(1 - phi^mal.m_MAL[1]) * inv(log10(1 - phi))
        sig = σ_fluid * phi^mal.m_MAL[1] + σ_matrix * (1 - phi)^p
    end

    return sig
end

function mix_models(σs, ϕ, ::single_phase)
    @assert length(σs) == 1
    return first(σs)
end

function mix_models(σs, ϕ, ::HS_plus_multi_phase)
    σ_min = minimum(σs)

    σ = inv(sum([ϕ[i] * inv(σs[i] + 2σ_min) for i in eachindex(σs)])) - 2σ_min

    return σ
end

function mix_models(σs, ϕ, ::HS_minus_multi_phase)
    σ_max = maximum(σs)

    σ = inv(sum([ϕ[i] * inv(σs[i] + 2σ_max) for i in eachindex(σs)])) - 2σ_max

    return σ
end

function mix_models(σs, ϕ, m::GAL)
    σ = sum([σs[i] * (ϕ[i]^m.m_GAL[i]) for i in eachindex(m.m_GAL)])
    return σ
end

function SubsurfaceCore.from_nt(m::T,
        ps::NamedTuple) where {T <: Union{
        HS1962_plus, HS1962_minus, HS_plus_multi_phase, HS_minus_multi_phase}}
    return m()
end

function SubsurfaceCore.from_nt(m::Type{T}, ps::NamedTuple) where {T <: MAL}
    return m(ps.m_MAL)
end

function SubsurfaceCore.from_nt(m::Type{T}, ps::NamedTuple) where {T <: GAL}
    return m(ps.m_GAL)
end

function sample_type(m::Type{T}) where {T <: phase_mixing}
    T
end
