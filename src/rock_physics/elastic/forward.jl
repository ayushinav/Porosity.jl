
function forward(m::anharmonic, p, params=default_params_anharmonic)
    @unpack T_K_ref, P_Pa_ref, Gu_0_ol, dG_dT, dG_dP, ν, Gu_0_crust, dG_dT_crust, dG_dP_crust, Gu_tp_fn, Ku_tp_fn = params

    Gu₀, dG_dT₀, dG_dP₀ = @. calc_Gu₀(
        Gu_0_ol, dG_dT, dG_dP, Gu_0_crust, dG_dT_crust, dG_dP_crust) #since χ is 1., we are always using ol

    ΔT = @. m.T - T_K_ref # K
    ΔP = @. m.P * 1.0f9 - P_Pa_ref # Pa
    Gu_TP = @. Gu_tp_fn(m.T, m.P, m.ρ)
    Ku_TP = @. Ku_tp_fn(m.T, m.P, m.ρ)
    Gu_tp = @. calc_Gu(Gu₀, ΔT, ΔP, dG_dT₀, dG_dP₀, Gu_TP)
    Ku_tp = @. calc_Ku(Gu_tp, ν, Ku_TP)

    Vp = @. calc_Vp(Ku_tp, Gu_tp, m.ρ)
    Vs = @. calc_Vs(Gu_tp, m.ρ)

    return RockphyElastic(Gu_tp, Ku_tp, Vp, Vs)
end

function forward(m::anharmonic_poro, p, params=default_params_anharmonic_poro)
    @unpack m_A, m_K, ν, p_anharmonic = params

    anh_p = forward(anharmonic(m.T, m.P, m.ρ), [], p_anharmonic)
    @unpack G, K = anh_p

    Γ_G = @. melt_shear_moduli(m.ϕ, m_A, ν)
    Gueff = @. Γ_G * G

    Γ_K = @. melt_bulk_moduli(m.ϕ, m_A, ν)
    K_sk = @. Γ_K * K

    nr = @. (1 - K_sk / K)^2
    dr = @. 1 - m.ϕ - K_sk / K + m.ϕ * K / m_K

    Kueff = @. K_sk + (nr / (dr + 1.0f-10)) * K
    μ_eff = @. G * Γ_G
    K_nr = @. (1 - Γ_K)^2
    K_dr = @. 1 - m.ϕ - Γ_K + m.ϕ * K / m_K
    K_eff = @. K * (Γ_K + K_nr / (K_dr + 1.0f-10))

    Vp = @. calc_Vp(K_eff, μ_eff, m.ρ)
    # @show μ_eff, K_eff
    Vs = @. calc_Vs(μ_eff, m.ρ)

    # return G, Γ_G, Γ_K, K_sk #RockphyElastic(Gueff, Kueff, Vp, Vs)
    return RockphyElastic(Gueff, Kueff, Vp, Vs)
end

function forward(m::SLB2005, p, params=(;))
    dV_P = @. 0.0380f0 * m.P
    dV_T = @. -0.000378f0 * (m.T - 300)

    Vs = @. 4.77 + dV_P + dV_T

    return RockphyElastic(zero(m.T .+ m.P), zero(m.T .+ m.P), zero(m.T .+ m.P), Vs * 1.0f3)
end

function forward(::Type{M}) where {M <: AbstractElasticModel}
    return RockphyElastic
end
