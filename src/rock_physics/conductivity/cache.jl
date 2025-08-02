#! format: off

# ========================================================================================================== 
# olivine

params_SEO3 = (
    S_bfe = 5.06f24, 
    H_bfe = 0.357f0,
    S_bmg = 4.58f26,
    H_bmg = 0.752f0,
    S_ufe = 12.2f-6,
    H_ufe = 1.05f0,
    S_umg = 2.72f-6,
    H_umg = 1.09f0,
)

params_UHO2014 = (
    # Ionic Vacancy
    H_v = 239f0, # kJ/mol, activation enthalpy
    S_v = 10f0 ^5.07f0, # S/m, pre-exponential conductivity

    # Polaron Hopping
    H_p = 144f0, # kJ/mol, activation enthalpy
    S_p = 10f0 ^2.34f0, # S/m, pre-exponential conductivity
    
    # Proton
    H_h = 89f0, # kJ/mol, activation enthalpy
    S_h = 10f0 ^-1.37f0, # S/m, pre-exponential conductivity
    a_h = 1.79f0, # (kJ/mol/wt)*(ppm ^1/3)
    r_h = 1f0, # unitless
)

params_Jones2012 = (
    S=10.0f0^(3.05f0),  # S/m, pre-exponential conductivity
    r = 0.86f0,  # unitless
    H = 0.91f0,  # eV, activation enthalpy
    a = 0.09f0,  # unitless
    params_SEO3 = params_SEO3
)

params_Poe2010 = (
    # Hydrous 100 axis
    S_H100 = 10f0^2.59f0,  # S/m, pre-exponential conductivity
    H_H100 = 1.26f0,  # eV, activation enthalpy
    a_H100 = 1.18f0,  # unitless
    r_H100 = 1f0,  # unitless
    Va_H100 = 0f0,  # cc/mol, activation volume

    # Hydrous 010 axis
    S_H010 = 10f0^3.46f0,  # S/m, pre-exponential conductivity
    H_H010 = 1.5f0,  # eV, activation enthalpy
    a_H010 = 1.43f0,  # unitless
    r_H010 = 1f0,  # unitless
    Va_H010 = 0f0,  # cc/mol, activation volume

    # Hydrous 001 axis
    S_H001 = 10f0^1.02f0,  # S/m, pre-exponential conductivity
    H_H001 = 0.812f0,  # eV, activation enthalpy
    a_H001 = 0.70f0,  # unitless
    r_H001 = 1f0,  # unitless
    Va_H001 = 0f0,  # cc/mol, activation volume

    # Anhydrous 100 axis
    S_A100 = 334f0,  # S/m, pre-exponential conductivity
    H_A100 = 1.46f0,  # eV, activation enthalpy
    Va_A100 = 0f0,  # cc/mol, activation volume

    # Anhydrous 010 axis
    S_A010 = 13.8f0,  # S/m, pre-exponential conductivity
    H_A010 = 1.12f0,  # eV, activation enthalpy
    Va_A010 = 0f0,  # cc/mol, activation volume

    # Anhydrous 001 axis
    S_A001 = 99f0,  # S/m, pre-exponential conductivity
    H_A001 = 1.29f0,  # eV, activation enthalpy
    Va_A001 = 0f0  # cc/mol, activation volume
)


params_Wang2006 = (
    # Hydrous
    S_H = 10f0^3.0f0,  # S/m, pre-exponential conductivity  
    H_H = 87f0,  # kJ/mol, activation enthalpy
    a_H = 0f0,  # unitless
    r_H = 0.62f0,  # unitless
    Va_H = 0f0,  # cc/mol, activation volume

    # Anhydrous
    S_A = 10f0^2.4f0,  # S/m, pre-exponential conductivity
    H_A = 154f0,  # kJ/mol, activation enthalpy
    Va_A = 0f0  # cc/mol, activation volume
)


const params_Yoshino2009 = (
    # Ionic conductivity parameters
    S_i = 10f0^4.73f0,  # S/m, pre-exponential conductivity
    H_i = 2.31f0,  # eV, activation enthalpy
    Va_i = 0f0,  # cc/mol, activation volume

    # Hopping Conduction parameters
    S_h = 10f0^2.98f0,  # S/m, pre-exponential conductivity
    H_h = 1.71f0,  # eV, activation enthalpy
    Va_h = 0f0,  # cc/mol, activation volume

    # Proton Conduction parameters
    S_p = 10f0^1.90f0,  # S/m, pre-exponential conductivity
    H_p = 0.92f0,  # eV, activation enthalpy
    a_p = 0.16f0,  # unitless
    r_p = 1f0,  # unitless
    Va_p = 0f0  # cc/mol, activation volume
)

# ========================================================================================================== 
# melts

params_Ni2011 = (
    T_corr = 1146.8f0, # K
    D = 0.006 # unitless, Partition coefficient {ol/melt}
)


params_Sifre2014 = (
    D_p = 0.007f0,  # unitless, D for peridotite/melt
    D_o = 0.002f0,  # unitless, D for ol/melt

    # Densities
    den_p = 3.3f0,  # g/cm^3, density peridotite
    den_h2o = 1.4f0,  # g/cm^3, density of water
    den_carb = 2.4f0,  # g/cm^3, density of molten carbonates [Liu and Lange, 2003]
    den_basalt = 2.8f0,  # g/cm^3, density of molten basalt [Lange and Carmichael, 1990]

    # H2O melt parameters
    a_h2o = 88774f0,
    b_h2o = 0.3880f0,
    c_h2o = 73029f0,
    d_h2o = 4.54f-5,
    e_h2o = 5.5607f0,

    # C2O melt parameters
    a_c2o = 789166f0,
    b_c2o = 0.1808f0,
    c_c2o = 32820f0,
    d_c2o = 5.50f-5,
    e_c2o = 5.7956f0
)
    
    
params_Gaillard2008 = (
    S = 3440, # S/m, pre-exponential conductivity
    H = 31.9 # KJ, activation enthalpy
)

# ========================================================================================================== 
# orthopyroxene

params_Dai_Karato2009 = (
    A = 10f0 ^ 2.4f0,
    Aw = 10f0 ^ 2.6f0,
    H = 147f0,
    Hw = 82f0,
    r = 62f-2
)

params_Zhang2012 = (
    S_pol = 3.99f0,
    H_pol = 1.88f0,
    S_hyd = 10f0 ^2.58f0,
    H_hyd = 84f-2,
    a = 8f-2,
)

# ========================================================================================================== 
# clinopyroxene

params_Yang2011 = (
    A = 10f0 ^ 2.16f0,
    Aw = 10f0 ^ 3.56f0,
    H = 1.06f0,
    Hw = 0.73f0,
    r= 113f-2
)

# olivine
const default_params_SEO3 = deepcopy(params_SEO3)
const default_params_UHO2014 = deepcopy(params_UHO2014)
const default_params_Jones2012 = deepcopy(params_Jones2012)
const default_params_Poe2010 = deepcopy(params_Poe2010)
const default_params_Wang2006 = deepcopy(params_Wang2006)
const default_params_Yoshino2009 = deepcopy(params_Yoshino2009)
# melts
const default_params_Ni2011 = deepcopy(params_Ni2011)
const default_params_Sifre2014 = deepcopy(params_Sifre2014)
const default_params_Gaillard2008 = deepcopy(params_Gaillard2008)
# orthopyroxene
const default_params_Dai_Karato2009 =deepcopy(params_Dai_Karato2009)
const default_params_Zhang2012 = deepcopy(params_Zhang2012)
#clinopyroxene
const default_params_Yang2011 =deepcopy(params_Yang2011)
