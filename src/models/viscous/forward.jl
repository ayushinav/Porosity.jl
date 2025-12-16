function SubsurfaceCore.forward(m::HZK2011, p, params=default_params_HZK2011)
    @unpack mechs, p_dep_calc, melt_enhancement = params

    P = @. p_dep_calc * m.P
    ϵ_rate = zero(m.T .+ m.P .+ m.dg .+ m.σ .+ m.ϕ) ## TODO
    x_ϕ_c_vec = get_melt_settings_for_x_ϕ_c(melt_enhancement)

    for mech in keys(mechs)
        sr = sr_flow_law_calculation(
            m.T, P, m.σ, m.dg, m.ϕ, 0, getfield(x_ϕ_c_vec, mech), getfield(mechs, mech))
        @. ϵ_rate += sr
    end

    η = @. m.σ * 1.0f3 * 1.0f6 / ϵ_rate

    return RockphyViscous(ϵ_rate, η)
end

function SubsurfaceCore.forward(m::HK2003, p, params=default_params_HK2003)
    @unpack mechs, ch2o_o, p_dep_calc, melt_enhancement = params

    P = @. p_dep_calc * m.P

    fH2O = @. calc_fH2O(m.Ch2o_ol, ch2o_o, P, m.T)

    ϵ_rate = zero(m.T .+ m.P .+ m.dg .+ m.σ .+ m.ϕ .+ m.Ch2o_ol)

    x_ϕ_c_vec = get_melt_settings_for_x_ϕ_c(melt_enhancement)

    for mech in keys(mechs)
        sr = broadcast(
            (T, P, σ, d, ϕ, fH2O) -> sr_flow_law_calculation_HK2003(
                T, P * 1.0f9, σ, d, ϕ, fH2O, getfield(x_ϕ_c_vec, mech), mechs, mech),
            m.T,
            P,
            m.σ,
            m.dg,
            m.ϕ,
            fH2O)
        @. ϵ_rate += sr
    end

    η = @. m.σ * 1.0f3 * 1.0f6 / ϵ_rate

    return RockphyViscous(ϵ_rate, η)
end

function SubsurfaceCore.forward(m::xfit_premelt, p, params=default_params_xfit_premelt)
    @unpack α, T_η, γ, B, Tr, Pr, η_r, H, V, M, dg_r = params

    Tprime = @. m.T / m.T_solidus
    A_n = @. calc_An(Tprime, m.ϕ, α, T_η, γ, B)

    η_meltfree = @. η_r *
                    (m.dg / dg_r)^M *
                    exp(V / (gas_R * 1.0f3) * (m.P / m.T - Pr / Tr) * 1.0f9 +
                        H / (gas_R * 1.0f3) * (1 / m.T - 1 / Tr))

    η = @. A_n * η_meltfree

    return RockphyViscous(zero(η), η)
end
