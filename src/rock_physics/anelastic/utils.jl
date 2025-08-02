# andrade_psp

function calc_X̃(T, d, P, ϕ, params_anelastic)
    @unpack n, β, τ_MR, E, G_UR, TR, PR, dR, Vstar, M, melt_alpha, ϕ_c, elastic_type, params_elastic, melt_enhancement = params_anelastic

    X̃ = @. (d / dR)^(-M) * exp((-E / (MT.gas_R * 1.0f3)) * (inv(T) - inv(TR)) -
                Vstar / (MT.gas_R * 1.0f3) * (P / T - PR / TR) * 1.0f9)
    x_ϕ_c = getfield(get_melt_settings_for_x_ϕ_c(Val{melt_enhancement}()), :diff)
    X̃ = @. X̃ / x_ϕ_c
    @. X̃ *= get_melt_enhancement(ϕ, melt_alpha, x_ϕ_c, ϕ_c)

    return X̃
end

function forward_for_anelastic(m, ::Val{anharmonic}, params)
    return forward(anharmonic(m.T, m.P, m.ρ), [], params)
end

function forward_for_anelastic(m, ::Val{anharmonic_poro}, params)
    return forward(anharmonic_poro(m.T, m.P, m.ρ, m.ϕ), [], params)
end

# eburgers_psp

function add_melt_effects(ϕ, scale, m_α, ϕ_c, x_ϕ_c)
    return scale * x_ϕ_c / get_melt_enhancement(ϕ, m_α, x_ϕ_c, ϕ_c)
end

function get_η_diff(m, viscous_type::Val{HZK2011}, params_viscous)
    @unpack mechs, p_dep_calc, melt_enhancement = params_viscous

    P = @. p_dep_calc * m.P
    x_ϕ_c_vec = get_melt_settings_for_x_ϕ_c(Val{melt_enhancement}())
    fH2O = 0.0f0
    ϵ_rate_diff = sr_flow_law_calculation(m.T, P * 1.0f9, m.σ, m.dg, m.ϕ, 0,
        getfield(x_ϕ_c_vec, :diff), getfield(mechs, :diff))
    η_diff = @. m.σ * 1.0f9 / ϵ_rate_diff

    return η_diff
end

function get_η_diff(m, viscous_type::Val{HK2003}, params_viscous)
    @unpack mechs, p_dep_calc, ch2o_o, melt_enhancement = params_viscous

    P = @. p_dep_calc * m.P
    x_ϕ_c_vec = get_melt_settings_for_x_ϕ_c(Val{melt_enhancement}())
    fH2O = @. calc_fH2O(m.Ch2o_ol, ch2o_o, P, m.T)
    ϵ_rate_diff = broadcast(
        (T, P, σ, d, ϕ, fH2O) -> sr_flow_law_calculation_HK2003(
            T, P * 1.0f9, σ, d, ϕ, fH2O, getfield(x_ϕ_c_vec, :diff), mechs, :diff),
        m.T,
        P,
        m.σ,
        m.dg,
        m.ϕ,
        fH2O)
    η_diff = @. m.σ * 1.0f9 / ϵ_rate_diff

    return η_diff
end

function get_η_diff(m, viscous_type::Val{xfit_premelt}, params_viscous)
    resp_xfit_premelt = forward(
        xfit_premelt(m.T, m.P, m.dg, m.σ, m.ϕ, m.T_solidus), [], params_viscous)

    return resp_xfit_premelt.η
    # requires T_solidus :)
end

function calc_maxwell_times(Gu, m::eburgers_psp, params_btype, JF10_visc,
        params_viscous, viscous_type, melt_enhancement)
    @unpack TR, PR, dR, E, Vstar, Tau_LR, Tau_HR, Tau_MR, Tau_PR, m_a, m_v, melt_alpha, ϕ_c = params_btype
    x_ϕ_c = getfield(get_melt_settings_for_x_ϕ_c(Val{melt_enhancement}()), :diff)

    if JF10_visc
        scale = @. (m.dg / dR)^m_v * exp(E / (MT.gas_R * 1.0f3) * (1 / m.T - 1 / TR) +
                       Vstar / (MT.gas_R * 1.0f3) * (m.P / m.T - PR / TR) * 1.0f9)
        new_scale = @. add_melt_effects(m.ϕ, scale, melt_alpha, ϕ_c, x_ϕ_c)
        τ_maxwell = @. Tau_MR * new_scale

    else
        # requires η_diff here
        η_diff = get_η_diff(m, Val{viscous_type}(), params_viscous)
        τ_maxwell = @. η_diff / Gu
    end

    LHP = @. (m.dg / dR)^m_a * exp(E / (MT.gas_R * 1.0f3) * (1 / m.T - 1 / TR) +
                 Vstar / (MT.gas_R * 1.0f3) * (m.P / m.T - PR / TR) * 1.0f9)
    scale_LHP = @. add_melt_effects(m.ϕ, LHP, melt_alpha, ϕ_c, x_ϕ_c)

    τ_L = @. Tau_LR * scale_LHP
    τ_H = @. Tau_HR * scale_LHP
    τ_P = @. Tau_PR * scale_LHP

    return τ_maxwell, τ_L, τ_H, τ_P
end

function integrate_fn(fn::T, a::T1, b::T2, N::Int, ::Val{:midpoint}) where {T, T1, T2}# Defining function for integrating using mid-point rule 
    dx = (b - a) / N
    mid_points = range(; start=a + dx / 2, stop=b - dx / 2, step=dx) # Mid-points of the intervals 
    f_vals = fn.(mid_points) # function value at the mid-points
    I = dx * sum(f_vals)
    return I
end

function integrate_fn(fn::T, a::T1, b::T2, N::Int, ::Val{:trapezoidal}) where {T, T1, T2} # Defining function for trapezoid rule 
    dx = (b - a) / N
    points = range(; start=a, stop=b, step=dx) # Edges of the interval
    f_vals = fn.(points) # function value at the edges
    I = 0
    for i in 1:N # Since the number of points will be N+1, going from 1 to N will include the first and last points only once
        I += (f_vals[i] + f_vals[i + 1]) / 2
    end
    I *= dx
    return I
end

function integrate_fn(fn::T, a::T1, b::T2, N::Int, ::Val{:simpson}) where {T, T1, T2} # Defining the function for Simpson's rule 
    dx = (b - a) / N
    # @show dx
    # points= collect(range(start= a, stop= b, step= dx/2)); # Getting the midpoints and the edges of all the intervals 
    points = collect(range(; start=a, stop=b, length=(2N + 1))) # Getting the midpoints and the edges of all the intervals 
    # points = [a, b]
    f_vals = @. fn(points)
    # f_vals= broadcast(fn, points);
    # @show length(f_vals), length(points)
    I = zero(a)
    for i in 0:(N - 1) # loop over all the intervals, the first interval is defined by i=0, second by i=1, and so on.
        j = 2i + 1 # Get the index position of the leftmost point of every interval in the array of all the points (defined as 'points')
        I += (f_vals[j] + 4f_vals[j + 1] + f_vals[j + 2]) / 6
        # I+= (fn(points[j])+ 4fn(points[j+1])+ fn(points[j+2]))/6 
    end
    I *= dx
    return I
end

function integrate_fn(fn, a::T1, b::T2, N::Int, ::Val{:quadgk}; kwargs...) where {T1, T2} # Defining the function for Simpson's rule dx= (b-a)/N;
    # @show a, b, N
    # @show quadgk(fn, a, b)
    return convert(typeof(a), first(quadgk(fn, a, b; kwargs...)))
end

function integrate_s(J_int_fn::F, ω::T, integration_params; kwargs...) where {F, T}
    @unpack l, h, N, rule = integration_params

    f(x) = J_int_fn(x, ω)
    # return integrate_fn(f, l, h, N, Val{:simpson}())
    return integrate_fn(f, l, h, N, Val{rule}(); kwargs...)
end

# function integrate(J_int_fn, ω::T, integration_params) where {T <: AbstractVector{<: Any}}
#     return [integrate(J_int_fn, iω, integration_params) for iω in vec(ω)]
# end

# function integrate(J_int_fn, ω::T, integration_params) where {T <: AbstractArray{<: Any}}
#     # @show ω
#     return [integrate(J_int_fn, iω, integration_params) for iω in vec(ω)]
# end

# xfit_premelt aka premelt_anelastic

function calc_Ap(Tn, ϕ, params)
    @unpack α_B, A_B, τ_pp, A_p_fac_1, A_p_fac_2, A_p_fac_3, σ_p_fac_1, σ_p_fac_2, σ_p_fac_3, A_p_Tn_pts, σ_p_Tn_pts, include_direct_melt_effect, β, β_B, poro_Λ = params

    β_p = (include_direct_melt_effect) ? β : 0.0f0

    A_p = 0.0f0
    if Tn >= A_p_Tn_pts[3]
        A_p = A_p_fac_3 + β_p * ϕ
    else
        if Tn >= A_p_Tn_pts[2]
            A_p = A_p_fac_3
        else
            if Tn >= A_p_Tn_pts[1]
                A_p = A_p_fac_1 + A_p_fac_2 * (Tn - A_p_Tn_pts[1])
            else
                A_p = A_p_fac_1
            end
        end
    end

    return A_p
end

function calc_σp(Tn, params)
    @unpack α_B, A_B, τ_pp, A_p_fac_1, A_p_fac_2, A_p_fac_3, σ_p_fac_1, σ_p_fac_2, σ_p_fac_3, A_p_Tn_pts, σ_p_Tn_pts, include_direct_melt_effect, β, β_B, poro_Λ = params

    σ_p = 0.0f0

    if Tn < σ_p_Tn_pts[1]
        σ_p = σ_p_fac_1
    elseif (Tn >= σ_p_Tn_pts[1]) && (Tn < σ_p_Tn_pts[2])
        σ_p = σ_p_fac_1 + σ_p_fac_2 * (Tn - σ_p_Tn_pts[1])
    else
        σ_p = σ_p_fac_3
    end

    return σ_p
end

# xfit_mxw

function get_fit_params(ps, fit)
    ps_fit = getfield(ps, fit)
    @unpack β2, τ_cutoff = ps_fit
    return β2, τ_cutoff
end

function xfit_mxw_func(τ, α_a, α_b, α_c, α2, β1, β2, α_τn, τ_cutoff)
    # @unpack fit, α_a, α_b, α_c, α_τn, α2, β1 = params
    # @unpack β2, τ_cutoff = getfield(params, fit)

    # α = α_a - α_b/(1 + α_c .* (τ_norm ^ α_τn))
    # β = ones(size(τ_norm)) .* β1
    # β[τ_norm < τ_cutoff] .= β2
    # α[τ_norm < τ_cutoff] .= α2

    # return β .* τ_norm.^α
    # @show τ, α_a, α_b, α_c, α2, β1, β2, α_τn, τ_cutoff

    β = (τ < τ_cutoff) ? β2 : β1
    α = (τ < τ_cutoff) ? α2 : α_a - α_b / (1.0f0 + α_c * (τ^α_τn))

    # β = β2
    # α = α2
    # @show α

    return β * τ^α
end
