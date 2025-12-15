#! format: off
const melt_settings_true = (diff=5.0f0, disl=1.0f0, gbs=2.5f0)
const melt_settings_false = (diff=1.0f0, disl=1.0f0, gbs=1.0f0)

# get_melt_settings_for_x_ϕ_c(::Val{true}) = melt_settings_true
# get_melt_settings_for_x_ϕ_c(::Val{false}) = melt_settings_false
get_melt_settings_for_x_ϕ_c(x::Bool) = x ? melt_settings_true : melt_settings_false

params_anharmonic = (
    Isaak1992=(T_K_ref=300, P_Pa_ref=1.0f5, Gu_0_ol=81, dG_dT=-13.6f6, dG_dP=1.8f0,
        ν=0.25f0, Gu_0_crust=40, dG_dT_crust=-3.0f6, dG_dP_crust=3.5f0,
        Gu_tp_fn=(T, P, ρ) -> -1.0f0, Ku_tp_fn=(T, P, ρ) -> -1.0f0),
    Cammarono2003=(
        T_K_ref=300, P_Pa_ref=1.0f5, Gu_0_ol=81, dG_dT=-14.0f6, dG_dP=1.4f0, ν=0.25f0,
        Gu_0_crust=40, dG_dT__crust=-3.0f6, dG_dP_crust=3.5f0, Gu_TP=-1.0f0, Ku_TP=-1.0f0))

params_anharmonic_poro = (
    m_A=1.6f0, m_K=30.0f9, ν=0.25f0, p_anharmonic=params_anharmonic.Isaak1992)

const default_params_anharmonic = deepcopy(params_anharmonic.Isaak1992)
const default_params_anharmonic_poro = deepcopy(params_anharmonic_poro)
