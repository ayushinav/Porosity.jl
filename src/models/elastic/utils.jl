# anharmonic

function calc_Gu₀(G1, dG_dT1, dG_dP1, G2, dG_dT2, dG_dP2; χ=1.0f0)
    Gu₀ = (G1 * χ + (1 - χ) * G2)
    dG_dT₀ = dG_dT1 * χ + (1 - χ) * dG_dT2
    dG_dP₀ = dG_dP1 * χ + (1 - χ) * dG_dP2

    return (Gu₀, dG_dT₀, dG_dP₀)
end

function calc_Gu(Gu_0, ΔT, ΔP, ∂G_∂T, ∂G_∂P, Gu=-1)
    if Gu < 0
        return Gu_0 * 1.0f9 + ΔT * ∂G_∂T + ΔP * ∂G_∂P
    else
        return Gu
    end
end

function calc_Ku(G, ν, Ku=-1)
    if Ku < 0
        return 2G / 3 * (1 + ν) / (1 - 2ν)
    else
        return Ku
    end
end

calc_Vp(K, G, ρ) = sqrt((K + 4G / 3) / ρ)
calc_Vs(G, ρ) = sqrt(G / ρ)

# anharmonic_poro

function melt_shear_moduli(ϕ::T1, A::T2, ν::T3) where {T1, T2, T3}
    ψ = 1 - A * sqrt(ϕ)

    # b = Float32[1.6122 0.13572 0.0
    #      4.5869 3.6086 0.0
    #      -7.5395 -4.8676 -4.3182]
    # ν_vec = ν .^ Float32[0, 1, 2]

    # b_vec = b * ν_vec # :)

    b_vec = Float32[1.64613, 5.4890504, -9.026288] # same as above, for type stability

    n_μ = b_vec[1] * ψ + b_vec[2] * (1 - ψ) + b_vec[3] * ψ * (1 - ψ)^2

    μ_sk_prime = (1 - (1 - ψ)^n_μ)
    Γ_G = (1 - ϕ) * μ_sk_prime

    return Γ_G
end

function melt_bulk_moduli(ϕ::T1, A::T2, ν::T3) where {T1, T2, T3}
    ψ = 1.0f0 - A * sqrt(ϕ)

    # a = Float32[1.8625 0.52594 -4.8397 0.0
    #               4.5001 -6.1551 -4.3634 0.0
    #               -5.6512 6.9159 29.595 -58.96]
    # ν_vec = ν .^ Float32[0, 1, 2, 3]

    # a_vec = a * ν_vec # :)

    a_vec = Float32[1.6915036, 2.6886127, -2.9937873] # same as above, for type stability

    n_k = a_vec[1] * ψ + a_vec[2] * (1 - ψ) + a_vec[3] * ψ * (1 - ψ)^(1.5f0)

    k_sk_prime = (1 - (1 - ψ)^n_k)
    Γ_K = (1 - ϕ) * k_sk_prime

    return Γ_K
end

function Vp_Vs_calc(ϕ, G, K, Γ_G, Γ_K, ρ, Km)
    μ_eff = G * Γ_G

    K_nr = (1 - Γ_K)^2
    K_dr = 1 - ϕ - Γ_K + ϕ * K / Km

    K_eff = K * (Γ_K + K_nr / (K_dr + 1.0f-10))

    Vp = calc_Vp(K_eff, μ_eff, ρ)
    Vs = calc_Vs(μ_eff, ρ)

    return Vp, Vs
end
