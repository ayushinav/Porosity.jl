"""
    tune_rp_modelType(fn_list, mtype)

## Arguments

  - `fn_list` : Vector of functions applied to the parameters in the same sequence
  - `mtype` : model type that will constructed using the given and estimated parameters

## Usage

```jldoctest
julia> T = (800:200:1200) .+ 273.0;

julia> P = 2.0;

julia> T_solidus = 1000.0 + 273;

julia> ps_nt = (; T=T, P=P, T_solidus=T_solidus);

julia> fn_list = [get_melt_fraction];

julia> m_type = two_phase_modelType(SEO3, Gaillard2008, HS1962_plus());

julia> m = tune_rp_modelType(fn_list, m_type)
Tuning rock physics with function list :
[Porosity.get_melt_fraction]
to obtain the rock physics model of type two_phase_modelType{SEO3, Gaillard2008, HS1962_plus}

julia> model = m(ps_nt)
Two phase composition using HS1962_plus()

*    m₁ (solid phase) : 
Model : SEO3
Temperature (K) : 1073.0:200.0:1473.0
ϕ : [1.0, 1.0, 0.4594594594593244]

*    m₂ (liquid phase) : 
Model : SEO3
Temperature (K) : 1073.0:200.0:1473.0
ϕ : [0.0, 0.0, 0.5405405405406756]
```
"""
mutable struct tune_rp_modelType{K, M} <: AbstractRockphyModel
    fn_list::K
    model::Type{M}
end

function tune_rp_modelType(fn_list, model)
    m_type = isa(model, Type) ? model : typeof(model)
    tune_rp_modelType(fn_list, m_type)
end

function (m::tune_rp_modelType{K, M})(ps_nt) where {K, M}
    for fn in m.fn_list
        nt = fn(ps_nt)
        ps_nt = (; ps_nt..., nt...)
    end
    return from_nt(M, ps_nt)
end

function (m::tune_rp_modelType{K, M})(ps_nt) where {K, M <: two_phase_modelType}
    for fn in m.fn_list
        nt = fn(ps_nt)
        ps_nt = (; ps_nt..., nt...)
    end
    return from_nt(m.model, ps_nt)
end

default_params(::Type{tune_rp_modelType{K, M}}) where {K, M} = default_params(M)
