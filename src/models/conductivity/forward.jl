# olivine

function SubsurfaceCore.forward(m::SEO3, p, params=default_params_SEO3)
    @unpack S_bfe, H_bfe, S_bmg, H_bmg, S_ufe, H_ufe, S_umg, H_umg = params
    fO₂ = fO2.(m.T)

    bfe = @. arrh_dry(S_bfe, H_bfe, boltz_k, m.T)
    bmg = @. arrh_dry(S_bmg, H_bmg, boltz_k, m.T)
    ufe = @. arrh_dry(S_ufe, H_ufe, boltz_k, m.T)
    umg = @. arrh_dry(S_umg, H_umg, boltz_k, m.T)

    concFe = @. bfe + 3.33f24 * exp(-0.02f0 * inv(boltz_k * m.T)) * fO₂^(1 / 6.0f0)
    concMg = @. bmg + 6.21f30 * exp(-1.83f0 * inv(boltz_k * m.T)) * fO₂^(1 / 6.0f0)
    σ = @. concFe * ufe * charge_e + 2.0f0 * concMg * umg * charge_e

    return RockphyCond(log10.(σ))
end

function SubsurfaceCore.forward(m::UHO2014, p, params=default_params_UHO2014)

    # Va = 0 so no dependence on pressure
    @unpack H_v, S_v, H_p, S_p, H_h, S_h, a_h, r_h = params

    σ_v = @. arrh_dry(S_v, H_v, gas_R, m.T)
    σ_p = @. arrh_dry(S_p, H_p, gas_R, m.T)
    σ_h = @. arrh_wet(S_h, H_h, gas_R, m.T, m.Ch2o_ol, a_h, r_h)

    σ = @. σ_v + σ_p + σ_h

    return RockphyCond(log10.(σ))
end

function SubsurfaceCore.forward(m::Jones2012, p, params=default_params_Jones2012)
    @unpack S, r, H, a, params_SEO3 = params

    @unpack S_bfe, H_bfe, S_bmg, H_bmg, S_ufe, H_ufe, S_umg, H_umg = params_SEO3

    # hydrous
    σ_H = @. arrh_wet(S, H, boltz_k, m.T, m.Ch2o_ol / 1.0f4, a, r)
    fO₂ = @. fO2(m.T)

    # anhydrous

    bfe = @. arrh_dry(S_bfe, H_bfe, boltz_k, m.T)
    bmg = @. arrh_dry(S_bmg, H_bmg, boltz_k, m.T)
    ufe = @. arrh_dry(S_ufe, H_ufe, boltz_k, m.T)
    umg = @. arrh_dry(S_umg, H_umg, boltz_k, m.T)
    concFe = @. bfe + 3.33f24 * exp(-0.02f0 * inv(boltz_k * m.T)) * fO₂^(1 / 6)
    concMg = @. bmg + 6.21f30 * exp(-1.83f0 * inv(boltz_k * m.T)) * fO₂^(1 / 6)
    σ_A = @. concFe * ufe * charge_e + 2.0f0 * concMg * umg * charge_e

    σ = @. σ_H + σ_A

    return RockphyCond(log10.(σ))
end

function SubsurfaceCore.forward(m::Poe2010, p, params=default_params_Poe2010)
    @unpack S_H100, H_H100, a_H100, r_H100, S_H010, H_H010, a_H010, r_H010, S_H001, H_H001,
    a_H001, r_H001, S_A100, H_A100, S_A010, H_A010, S_A001, H_A001 = params

    # Anhydrous
    σ_A100 = @. arrh_dry(S_A100, H_A100, boltz_k, m.T)
    σ_A010 = @. arrh_dry(S_A010, H_A010, boltz_k, m.T)
    σ_A001 = @. arrh_dry(S_A001, H_A001, boltz_k, m.T)
    σ_A = @. (σ_A001 * σ_A010 * σ_A100)

    # Hydrous
    σ_H100 = @. arrh_wet(S_H100, H_H100, boltz_k, m.T, m.Ch2o_ol * 1.0f-4, a_H100, r_H100)
    σ_H010 = @. arrh_wet(S_H010, H_H010, boltz_k, m.T, m.Ch2o_ol * 1.0f-4, a_H010, r_H010)
    σ_H001 = @. arrh_wet(S_H001, H_H001, boltz_k, m.T, m.Ch2o_ol * 1.0f-4, a_H001, r_H001)
    σ_H = @. cbrt(σ_H001 * σ_H010 * σ_H100)

    σ = @. σ_H + σ_A

    return RockphyCond(log10.(σ))
end

function SubsurfaceCore.forward(m::Yoshino2009, p, params=default_params_Yoshino2009)
    @unpack S_i, H_i, S_h, H_h, S_p, H_p, a_p, r_p = params

    # ionic conduction
    σ_i = @. arrh_dry(S_i, H_i, boltz_k, m.T)

    # polaron hopping
    σ_h = @. arrh_dry(S_h, H_h, boltz_k, m.T)

    # proton conduction
    σ_p = @. arrh_wet(S_p, H_p, boltz_k, m.T, m.Ch2o_ol * 1.0f-4, a_p, r_p)

    σ = @. σ_i + σ_h + σ_p

    return RockphyCond(log10.(σ))
end

function SubsurfaceCore.forward(m::Wang2006, p, params=default_params_Wang2006)
    @unpack S_H, H_H, a_H, r_H, S_A, H_A = params

    σ_A = @. arrh_dry(S_A, H_A, gas_R, m.T)
    σ_H = @. arrh_wet(S_H, H_H, gas_R, m.T, m.Ch2o_ol * 1.0f-4, a_H, r_H)

    σ = @. σ_H + σ_A

    return RockphyCond(log10.(σ))
end

function SubsurfaceCore.forward(m::const_matrix, p, params=default_params_const_matrix)
    RockphyCond(log10.(m.cond))
end

# ========================================================================================================== 
# melt

function SubsurfaceCore.forward(m::Ni2011, p, params=default_params_Ni2011)
    @unpack T_corr, D = params

    ls = @. 2.172f0 - (860.82f0 - 204.46f0 * sqrt(m.Ch2o_m / 1.0f4)) * inv(m.T - T_corr)
    σ = @. 10.0f0^ls

    return RockphyCond(log10.(σ))
end

function SubsurfaceCore.forward(m::Sifre2014, p, params=default_params_Sifre2014)
    @unpack a_h2o, b_h2o, c_h2o, d_h2o, e_h2o, a_c2o, b_c2o, c_c2o, d_c2o, e_c2o = params

    H_h2o = @. a_h2o * exp(-b_h2o * m.Ch2o_m * 1.0f-4) + c_h2o
    lS_h2o = @. d_h2o * H_h2o + e_h2o
    S_h2o = @. exp(lS_h2o)
    melt_h2o = @. S_h2o * exp(-H_h2o * inv(1.0f3 * gas_R * m.T))

    # Esig CO2 melt
    H_co2 = @. a_c2o * exp(-b_c2o * m.Cco2_m * 1.0f-4) + c_c2o
    lS_co2 = @. d_c2o * H_co2 + e_c2o
    S_co2 = @. exp(lS_co2)
    melt_co2 = @. S_co2 * exp(-H_co2 * inv(1.0f3 * gas_R * m.T))

    # summation of conduction mechanisms
    σ = @. melt_co2 + melt_h2o

    return RockphyCond(log10.(σ))
end

function SubsurfaceCore.forward(m::Gaillard2008, p, params=default_params_Gaillard2008)
    @unpack S, H = params
    σ = @. arrh_dry(S, H, gas_R, m.T)
    return RockphyCond(log10.(σ))
end

# ========================================================================================================== 
# orthopyroxene

function SubsurfaceCore.forward(m::Dai_Karato2009, p, params=default_params_Dai_Karato2009)
    @unpack A, Aw, H, Hw, r = params

    σ_dry = @. arrh_dry(A, H, gas_R, (m.T - 273))
    σ_wet = @. arrh_wet(Aw, Hw, gas_R, (m.T - 273), m.Ch2o_opx * 1.0f-4, 0.0f0, r)

    σ = @. σ_dry + σ_wet

    return RockphyCond(log10.(σ))
end

function SubsurfaceCore.forward(m::Zhang2012, p, params=default_params_Zhang2012)
    @unpack S_pol, H_pol, S_hyd, H_hyd, a = params

    σ_pol = @. arrh_dry(S_pol, H_pol, boltz_k, m.T)
    σ_hyd = @. arrh_wet(S_hyd, H_hyd, boltz_k, m.T, m.Ch2o_opx * 1.0f-4, a, 1.0f0)

    σ = @. σ_pol + σ_hyd

    return RockphyCond(log10.(σ))
end

# ========================================================================================================== 
# clinopyroxene

function SubsurfaceCore.forward(m::Yang2011, p, params=default_params_Yang2011)
    @unpack A, Aw, H, Hw, r = params

    σ_dry = @. arrh_dry(A, H, boltz_k, m.T)
    σ_wet = @. arrh_wet(Aw, Hw, boltz_k, m.T, m.Ch2o_cpx * 1.0f-4, 0.0f0, r)

    σ = @. σ_dry + σ_wet

    return RockphyCond(log10.(σ))
end
