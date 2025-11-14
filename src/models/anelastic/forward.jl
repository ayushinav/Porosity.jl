function SubsurfaceCore.forward(m::andrade_psp, p, params=default_params_andrade_psp)
    @unpack n, β, τ_MR, E, G_UR, TR, PR, dR, Vstar, M, melt_alpha, ϕ_c, elastic_type, params_elastic, melt_enhancement = params

    resp_elastic = forward_for_anelastic(m, Val{elastic_type}(), params_elastic)
    @unpack G, K, Vp, Vs = resp_elastic

    Ju = @. inv(G)
    ω = Float32(2π) .* m.f

    X̃ = calc_X̃(m.T, m.dg, m.P, m.ϕ, params)

    param1 = @. β * gamma(1.0f0 + n) * cos(π / 2.0f0 * n)
    param2 = @. β * gamma(1.0f0 + n) * sin(π / 2.0f0 * n)

    ωX = @. ω ./ X̃

    J1 = @. Ju * (1 + param1 * inv(ωX)^n)
    J2 = @. Ju * (param2 * inv(ωX)^n + inv.(τ_MR * ωX))
    Qinv = @. J2 * inv(J1)
    Ma = @. sqrt.(inv(J1^2 + J2^2))
    Va = @. sqrt.(Ma / m.ρ)

    d_ = length(size(Va))
    Vave = dropdims(sum(Va; dims=d_); dims=d_) * Float32(inv(length(m.f)))

    return RockphyAnelastic(J1, J2, Qinv, Ma, Va, Vave)
end

function SubsurfaceCore.forward(m::eburgers_psp, p, params=default_params_eburgers_psp)
    @unpack integration_params, elastic_type, params_elastic, params_btype, viscous_type, params_viscous, JF10_visc, melt_enhancement = params
    @unpack alf, DeltaB, DeltaP, sig = params_btype

    ω = 2.0f0π .* m.f

    resp_elastic = resp_elastic = forward_for_anelastic(
        m, Val{elastic_type}(), params_elastic)
    @unpack G, K, Vp, Vs = resp_elastic

    Ju = @. inv(G)

    τ_maxwell, τ_L, τ_H, τ_P = calc_maxwell_times(
        G, m, params_btype, JF10_visc, params_viscous,
        Symbol(viscous_type), melt_enhancement)

    J1_int_fn(x, ω) = x^(alf - 1) / (1 + (ω * x)^2)
    J2_int_fn(x, ω) = x^alf / (1 + (ω * x)^2)

    int1 = broadcast(
        (l, h, omega) -> integrate_s(J1_int_fn,
            omega,
            (l=l, h=h, N=integration_params.τ_integration_points,
                rule=integration_params.integration_method)),
        τ_L,
        τ_H,
        ω)
    int2 = broadcast(
        (l, h, omega) -> integrate_s(J2_int_fn,
            omega,
            (l=l, h=h, N=integration_params.τ_integration_points,
                rule=integration_params.integration_method)),
        τ_L,
        τ_H,
        ω)

    τ_fac = @. alf * DeltaB / (τ_H^alf - τ_L^alf)

    J1 = @. 1 + τ_fac * int1
    J2 = @. ω * τ_fac * int2 + inv(ω * τ_maxwell)

    if DeltaP > 0 # TODO :check for bugs (?)
        function J1_int_fn2(x, ω, tau_p)
            inv(x) * (exp(-0.5f0 * log(x / tau_p) * inv(sig))^2) * inv(1 + (ω * x)^2)
        end
        int11 = broadcast(
            (omega, tau_p) -> integrate_s((x, omega) -> J1_int_fn2(x, omega, tau_p), omega,
                (l=eps(typeof(omega))^40, h=Inf, N=1, rule=:quadgk);
                rtol=1.0f16 * eps(eltype(ω))),
            ω,
            τ_P) # TODO : check this case

        @. J1 = J1 + DeltaP * int11 * inv(sig * sqrt(2.0f0π))

        function J2_int_fn2(x, ω, tau_p)
            (exp(-0.5f0 * log(x / tau_p) * inv(sig))^2) * inv(1 + (ω * x)^2)
        end
        int22 = broadcast(
            (omega, tau_p) -> integrate_s((x, omega) -> J2_int_fn2(x, omega, tau_p),
                omega, (l=0.0f0, h=Inf, N=1, rule=:quadgk)),
            ω,
            τ_P) # TODO : check this case

        @. J2 = J2 + DeltaP * ω * int22 * inv(sig * sqrt(2.0f0π))
    end

    @. J1 = Ju * J1
    @. J2 = Ju * J2

    Qinv = @. J2 * inv(J1)
    Ma = sqrt.(inv.(J1 .^ 2 + J2 .^ 2))
    Va = sqrt.(Ma ./ m.ρ)

    d_ = length(size(Va))
    Vave = dropdims(sum(Va; dims=d_); dims=d_) * Float32(inv(length(m.f)))

    return RockphyAnelastic(J1, J2, Qinv, Ma, Va, Vave)
end

function SubsurfaceCore.forward(
        m::premelt_anelastic, p, params=default_params_premelt_anelastic)
    @unpack params_xfit, elastic_type, elastic_params, viscous_params = params
    @unpack include_direct_melt_effect, β_B, poro_Λ, α_B, A_B, τ_pp = params_xfit

    resp_elastic = forward_for_anelastic(m, Val{elastic_type}(), elastic_params)
    @unpack G, K, Vp, Vs = resp_elastic

    Ju = @. inv(G)

    Tn = @. m.T / m.T_solidus

    viscous_type = xfit_premelt
    η = get_η_diff(m, Val{viscous_type}(), viscous_params)

    τ_m = @. η * Ju

    Apfn(Tn, ϕ) = calc_Ap(Tn, ϕ, params_xfit)
    A_p = broadcast(Apfn, Tn, m.ϕ)

    σpfn(Tn) = calc_σp(Tn, params_xfit)
    σ_p = broadcast(σpfn, Tn)

    β_B = @. (include_direct_melt_effect) ? β_B * m.ϕ : 0.0f0
    poro_elastic_factor = @. (include_direct_melt_effect) ? poro_Λ * m.ϕ : 0.0f0
    k_temp = @. A_B + β_B
    pifac = sqrt(π / 2.0f0)

    T = @. inv(m.f)
    p_p = @. T / (2.0f0π * τ_m)
    ABppa = @. k_temp * (p_p)^(α_B)
    lntauapp = @. log(τ_pp / p_p)

    J1 = @. Ju * (1 +
             poro_elastic_factor +
             ABppa / α_B +
             pifac * A_p * σ_p * erfc(lntauapp / (sqrt(2.0f0) * σ_p)))
    J2 = @. Ju * π / 2.0f0 * (ABppa + A_p * exp(-((lntauapp .^ 2) / (2 * σ_p^2))))

    Qinv = @. J2 * inv(J1)
    Ma = sqrt.(inv.(J1 .^ 2 + J2 .^ 2))
    Va = sqrt.(Ma ./ m.ρ)

    d_ = length(size(Va))
    Vave = dropdims(sum(Va; dims=d_); dims=d_) * Float32(inv(length(m.f)))

    return RockphyAnelastic(J1, J2, Qinv, Ma, Va, Vave)
end

function SubsurfaceCore.forward(m::xfit_mxw, p, params=default_params_xfit_mxw)
    @unpack α_a, α_b, α_c, α_τn, α2, β1, β2, τ_cutoff, melt_alpha, ϕ_c, elastic_type, elastic_params, viscous_type, viscous_params = params

    resp_elastic = forward_for_anelastic(m, Val{elastic_type}(), elastic_params)
    @unpack G, K, Vp, Vs = resp_elastic

    Ju = @. inv(G)

    ω = @. 2.0f0π * m.f
    τ = @. inv(ω)

    η_diff = get_η_diff(m, Val{viscous_type}(), viscous_params)
    τ_maxwell = @. η_diff / G

    τ_norm = @. τ / τ_maxwell
    f_norm = @. τ_maxwell * m.f

    τ_norm_f = @. inv(2.0f0π * f_norm)

    J_int_fn(x, _) = inv(x) * xfit_mxw_func(x, α_a, α_b, α_c, α2, β1, β2, α_τn, τ_cutoff)

    int1 = broadcast(
        tau_norm_f -> integrate_s(
            J_int_fn, 0.0f0, (l=10.0f0^(-30.0f0), h=tau_norm_f, N=1, rule=:quadgk)),
        τ_norm_f) # TODO : check this case

    int2 = broadcast(J_int_fn, τ_norm_f, 0.0f0)

    J1 = @. Ju * (1.0f0 + int1)
    J2 = @. Ju * (π / 2.0f0 * int2 * (τ_norm_f) + τ_norm)

    Qinv = @. J2 * inv(J1)
    Ma = sqrt.(inv.(J1 .^ 2 + J2 .^ 2))
    Va = sqrt.(Ma ./ m.ρ)

    d_ = length(size(Va))
    Vave = dropdims(sum(Va; dims=d_); dims=d_) * Float32(inv(length(m.f)))

    return RockphyAnelastic(J1, J2, Qinv, Ma, Va, Vave)
end

function SubsurfaceCore.forward(
        m::andrade_analytical, p, params=default_params_andrade_analytical)
    @unpack α, β, η_ss, viscosity_method, viscosity_mech, elastic_type, elastic_params, viscous_type, viscous_params = params

    resp_elastic = forward_for_anelastic(m, Val{elastic_type}(), elastic_params)
    @unpack G, K, Vp, Vs = resp_elastic

    Ju = @. inv(G)
    ω = @. 2.0f0π * m.f

    if viscosity_method
        η = get_η_diff(m, Val{viscous_type}(), viscous_params) # CHANGE HERE for any mechanism
    else
        η = η_ss
    end

    τ_maxwell = @. η / G

    MJ_real = @. 1 + β * gamma(1 + α) * cos(α * π / 2.0f0) * inv.(α * ω)
    MJ_imag = @. inv(ω * τ_maxwell) + β * gamma(1 + α) * sin(α * π / 2.0f0) * inv(ω^α)

    J1 = @. Ju * MJ_real
    J2 = @. Ju * MJ_imag

    J1J2_fac = @. 0.5f0 * (1 + sqrt(1 + (J2 / J1)^2))
    Qinv = J2 ./ J1 .* J1J2_fac

    Ma = sqrt.(inv.(J1 .^ 2 + J2 .^ 2))
    Va = sqrt.(Ma ./ m.ρ)

    d_ = length(size(Va))
    Vave = dropdims(sum(Va; dims=d_); dims=d_) * Float32(inv(length(m.f)))

    return RockphyAnelastic(J1, J2, Qinv, Ma, Va, Vave)
end
