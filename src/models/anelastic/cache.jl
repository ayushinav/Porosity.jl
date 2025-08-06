#! format: off
params_andrade_psp = (
    n = inv(3f0),
    β = 2f-2,
    τ_MR = 10f0 ^ (5.3f0),
    E = 303f3,
    G_UR = 62.2f0, # GPa,
    TR = 1173f0,
    PR = 2f-1,
    dR = 3.1f0,
    Vstar = 1f-5,
    M = 1f0,
    melt_alpha = 25,
    ϕ_c = 1f-5,
    # x_ϕ_c = 5,
    elastic_type = anharmonic,
    params_elastic = params_anharmonic.Isaak1992,
    melt_enhancement = false
)


params_JF10 = (
    bg_only = (
        dR = 13.4f0,  # ref grain size in microns
        G_UR = 62.5f0,  # GPa, unrel. G, reference val.
        E = 303000f0,  # J/mol
        m_a = 1.19f0, # grain size exponent for tau_i,i in (L,H,P)
        alf = 0.257f0, # high temp background tau exponent
        DeltaB = 1.13f0, # relaxation strength
        Tau_LR = 1f-3,  # Relaxation time lower limit reference
        Tau_HR = 1f7,  # Relaxation time higher limit reference
        Tau_MR = 10^6.95f0, # Reference Maxwell relaxation time
        DeltaP = 0f0, # no peak, set to 0
        sig = 0f0, # no peak, set to 0
        Tau_PR = 0f0, # no peak, set to 0
        TR = 1173f0,
        PR = 0.2f0,
        Vstar = 1f-5,
        m_v = 3f0, 
        melt_alpha = 25f0,
        ϕ_c = 1f-5,
        # x_ϕ_c = 5,
    ),
    bg_peak = (
        DeltaP = 0.057f0,  # relaxation strength of peak
        sig = 4f0,  # sigma, peak breadth
        Tau_PR = 10^-3.4f0,  # center maxwell time
        dR = 13.4f0,  # ref grain size in microns
        G_UR = 66.5f0,  # GPa, unrel. G, reference val.
        E = 360000f0,  # J/mol
        m_a = 1.31f0,  # grain size exponent for tau_i, i in (L,H,P)
        alf = 0.274f0,  # high temp background tau exponent
        DeltaB = 1.13f0,  # relaxation strength of background
        Tau_LR = 1f-3,  # Relaxation time lower limit reference
        Tau_HR = 1f7,  # Relaxation time higher limit reference
        Tau_MR = 10^7.48f0,  # Reference Maxwell relaxation time
        TR = 1173f0,
        PR = 0.2f0,
        Vstar = 1f-5,
        m_v = 3f0,
        melt_alpha = 25f0,
        ϕ_c = 1f-5,
        # x_ϕ_c = 5,
    ),
    s6585_bg_only = (
        dR = 3.1f0,  # ref grain size in microns
        G_UR = 62.0f0,  # GPa, unrel. G, reference val.
        E = 303000f0,  # J/mol
        m_a = 1.19f0,  # grain size exponent for tau_i, i in (L,H,P)
        alf = 0.33f0,  # high temp background tau exponent
        DeltaB = 1.4f0,  # relaxation strength
        Tau_LR = 1f-2,  # Relaxation time lower limit reference
        Tau_HR = 1f6,  # Relaxation time higher limit reference
        Tau_MR = 10^5.2f0,  # Reference Maxwell relaxation time
        DeltaP = 0f0,  # no peak, set to 0
        sig = 0f0,  # no peak, set to 0
        Tau_PR = 0f0,  # no peak, set to 0
        TR = 1173f0,
        PR = 0.2f0,
        Vstar = 1f-5,
        m_v = 3f0,
        melt_alpha = 25f0,
        ϕ_c = 1f-5,
        # x_ϕ_c = 5,
    ),
    s6585_bg_peak = (
        DeltaP = 0.07f0,  # relaxation strength of peak
        sig = 4f0,  # sigma, peak breadth
        Tau_PR = 10^-2.9f0,  # center maxwell time
        dR = 3.1f0,  # ref grain size in microns
        G_UR = 66.5f0,  # GPa, unrel. G, reference val.
        E = 327000f0,  # J/mol
        m_a = 1.19f0,  # grain size exponent for tau_i, i in (L,H,P)
        alf = 0.33f0,  # high temp background tau exponent
        DeltaB = 1.4f0,  # relaxation strength
        Tau_LR = 1f-2,  # Relaxation time lower limit reference
        Tau_HR = 1f6,  # Relaxation time higher limit reference
        Tau_MR = 10^5.4f0,  # Reference Maxwell relaxation time
        TR = 1173f0,
        PR = 0.2f0,
        Vstar = 1f-5,
        m_v = 3f0,
        melt_alpha = 25f0,
        ϕ_c = 1f-5,
        # x_ϕ_c = 5,
    )
)

params_eburgers_psp = (
    integration_params = (
        nτ = 3000,
        integration_method = :quadgk,
        τ_integration_points = 5000
    ),

    params_btype = params_JF10.bg_only,
    elastic_type = anharmonic,
    params_elastic = params_anharmonic.Isaak1992,
    viscous_type = xfit_premelt,
    params_viscous = params_xfit_premelt,
    melt_enhancement = false,
    JF10_visc = true

)

params_premelt_anelastic = (
    params_xfit = ( 
        α_B = 0.38f0,  # high temp background exponent
        A_B = 0.664f0,  # high temp background dissipation strength

        # Prf-melting dissipation peak settings
        τ_pp = 6f-5,  # peak center, table 4 of YT16, paragraph before eq 10
        A_p_fac_1 = 0.01f0,
        A_p_fac_2 = 0.4f0,
        A_p_fac_3 = 0.03f0,
        σ_p_fac_1 = 4f0,
        σ_p_fac_2 = 37.5f0,
        σ_p_fac_3 = 7f0,
        A_p_Tn_pts = [0.91f0, 0.96f0, 1f0],  # Tn cutoff points
        σ_p_Tn_pts = [0.92f0, 1f0],  # Tn cutoff points

        # Melt effects
        include_direct_melt_effect = false,  # set to true to include YT2024 melt effect
        β = 1.38f0,  # determined in YT2024, named Beta_P in YT2024 eq 5
        β_B = 6.94f0,  # YT2024 only
        poro_Λ = 4.0f0, # Table 6 YT2024
    ),

    elastic_type = anharmonic,
    elastic_params = params_anharmonic.Isaak1992,

    viscous_type = xfit_premelt,
    viscous_params = params_xfit_premelt
)

params_xfit_mxw = (
    fit1 = (
        β2 = 1853f0,
        τ_cutoff = 1f-11,
        α2 = 0.5f0,
        β1 = 0.32f0,
        α_a = 0.39f0,
        α_b = 0.28f0,
        α_c = 2.6f0,
        α_τn = 1f-1,
        melt_alpha = 25f0,
        ϕ_c = 1f-5,

        elastic_type = anharmonic,
        elastic_params = params_anharmonic.Isaak1992,

        viscous_type = xfit_premelt,
        viscous_params = params_xfit_premelt
    ),

    fit2 = (
        β2 = 8.476f0,
        τ_cutoff = 5f-6,
        α2 = 0.5f0,
        β1 = 0.32f0,
        α_a = 0.39f0,
        α_b = 0.28f0,
        α_c = 2.6f0,
        α_τn = 1f-1,
        melt_alpha = 25f0,
        ϕ_c = 1f-5,

        elastic_type = anharmonic,
        elastic_params = params_anharmonic.Isaak1992,

        viscous_type = xfit_premelt,
        viscous_params = params_xfit_premelt
    ),
)

params_andrade_analytical = (
    α = 1/3f0,
    β =1f-4,
    η_ss = 1f23,
    viscosity_method = true,
    viscosity_mech = :diff,

    elastic_type = anharmonic,
    elastic_params = params_anharmonic.Isaak1992,

    viscous_type = HK2003,
    viscous_params = params_HK2003
)

const default_params_andrade_psp = deepcopy(params_andrade_psp)
const default_params_eburgers_psp = deepcopy(params_eburgers_psp)
const default_params_premelt_anelastic = deepcopy(params_premelt_anelastic)
const default_params_xfit_mxw = deepcopy(params_xfit_mxw.fit1)
const default_params_andrade_analytical = deepcopy(params_andrade_analytical)
