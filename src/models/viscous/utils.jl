# HZK2011

function get_melt_enhancement(phi, α, x_ϕ_c, ϕ_c)
    a = log(x_ϕ_c)
    ratefac = inv(ϕ_c)
    step = a * erf(phi * ratefac)
    slope = α * phi
    ln_SR_phi_enh = slope + step
    SR_phi_enh = exp(ln_SR_phi_enh)

    return SR_phi_enh
end

function sr_flow_law_calculation(T, P, σ, d, ϕ, fH2O, x_ϕ_c, params)
    # @show params
    @unpack A, Q, V, p, n, alf, r, ϕ_c = params

    # (x_phi_c == 0) && (x_ϕ_c = 1f0)

    sr = @. inv(x_ϕ_c) *
            A *
            ((σ * 1.0f3)^n) *
            (d^(-p)) *
            exp(-(Q + P * 1.0f9 * V) / (gas_R * T * 1.0f3)) *
            (fH2O^r)
    enhance = @. get_melt_enhancement(ϕ, alf, x_ϕ_c, ϕ_c)

    return @. sr * enhance
end

function sr_flow_law_calculation_HK2003(T, P, σ, d, ϕ, fH2O, x_ϕ_c, mechs, mech)
    @unpack A, Q, V, p, n, alf, r, ϕ_c = HK2003_mech(T, fH2O, mechs, mech)

    sr = @. inv(x_ϕ_c) *
            A *
            ((σ * 1.0f3)^n) *
            (d^(-p)) *
            exp(-(Q + P * V) / (gas_R * T * 1.0f3)) *
            (fH2O^r)
    enhance = @. get_melt_enhancement(ϕ, alf, x_ϕ_c, ϕ_c)

    return @. sr * enhance
end

# HK2003

function calc_fH2O(H2O_ppm, H2O_o, P, T)
    E = 40.0f3
    V = 10.0f-6
    A_o = 26

    return (H2O_ppm >= H2O_o) * H2O_ppm / A_o * exp((E + P * V) / (gas_R * 1.0f3 * T))
end

# function HK2003_mech2!(pvec, T, fH2O, mechs, mech)
#     if mech == :diff
#         if fH2O > 0
#             pvec .=  getfield(getfield(mechs, mech), :wet) |> values
#         else
#             pvec .=  getfield(getfield(mechs, mech), :dry) |> values
#         end
#     elseif mech == :disl
#         if fH2O > 0
#             pvec .=  getfield(getfield(mechs, mech), :wet) |> values
#         else
#             pvec .=  getfield(getfield(mechs, mech), :dry) |> values
#         end
#     elseif mech == :gbs
#         if T >= (1250 + 273)
#             pvec .=  getfield(getfield(mechs, mech), :gt1250) |> values
#         else
#             pvec .=  getfield(getfield(mechs, mech), :lt1250) |> values
#         end
#     end

#     # x_ϕ_c_vec = get_melt_settings_for_x_ϕ_c(Val{melt_enhancement}())

#     return pvec #(ps..., x_ϕ_c = getfield(x_ϕ_c_vec, mech))
# end

function HK2003_mech(T, fH2O, mechs, mech)
    if mech == :diff
        if fH2O > 0
            ps = getfield(getfield(mechs, mech), :wet)
        else
            ps = getfield(getfield(mechs, mech), :dry)
        end
    elseif mech == :disl
        if fH2O > 0
            ps = getfield(getfield(mechs, mech), :wet)
        else
            ps = getfield(getfield(mechs, mech), :dry)
        end
    elseif mech == :gbs
        if T >= (1250 + 273)
            ps = getfield(getfield(mechs, mech), :gt1250)
        else
            ps = getfield(getfield(mechs, mech), :lt1250)
        end
    end

    return ps
end

# xfit_premelt

function calc_An(Tprime, ϕ, α, T_η, γ, B)
    if Tprime < T_η
        return 1.0f0
    else
        if (Tprime < 1)
            return exp(-log(γ) * (Tprime - T_η) / (Tprime * (1.0f0 - T_η)))
        else
            # @show α, ϕ, γ, B
            return exp(-α * ϕ) * inv(γ * B)
        end
    end
end
