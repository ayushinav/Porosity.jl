#! format: off
params_HZK2011 = (
    mechs = (
        diff=(
            A=10.0f0^(76.0f-1),  # Preexponential for coble diffusion creep
            Q=375.0f3,   # Activation energy for coble diffusion creep
            V=10.0f-6,   # Activation volume for coble diffusion creep
            p=3,       # Grain size exponent
            alf=25,    # Melt factor
            r=0,       # Water fugacity exponent
            n=1,       # Stress exponent
            ϕ_c=1.0f-5,
            # x_ϕ_c=5f0
        ), 
        
        disl=(
            A=1.1f5,   # Preexponential
            Q=530.0f3,   # Activation energy
            V=15.0f-6,   # Activation volume
            n=35.0f-1,     # Stress exponent
            p=0,       # Grain size exponent
            alf=30,    # Melt factor
            r=0,       # Water fugacity exponent
            ϕ_c=1.0f-5,
            # x_ϕ_c=1f0
        ), 
        
        gbs=(
            A=10.0f0^4.8f0,  # Preexponential for GBS-disl creep
            Q=445.0f3,   # Activation energy for GBS-disl creep
            V=15.0f-6,   # Activation volume
            p=73.0f-2,    # Grain size exponent
            n=29.0f-1,     # Stress exponent
            alf=35,    # Melt factor
            r=0,       # Water fugacity exponent
            ϕ_c=1.0f-5,
            # x_ϕ_c=25.0f-1
        ),
        ),
    p_dep_calc = true,
    melt_enhancement = false
    
)

params_HK2003 = (
    mechs = (
        diff = (
            dry = (
                A = 1.5f9,   # Preexponential for coble diffusion creep
                Q = 375.0f3, # Activation energy for coble diffusion creep
                V = 10.0f-6, # Activation volume for coble diffusion creep
                p = 3f0,       # Grain size exponent
                n = 1f0,       # Stress exponent
                r = 0f0,       # Water fugacity exponent
                alf = 25f0,    # Melt factor
                ϕ_c = 1.0f-5,
                # x_ϕ_c = 5f0
            ),
            wet = (
                A = 2.5f7,   # Preexponential for coble diffusion creep
                Q = 375.0f3, # Activation energy for coble diffusion creep
                V = 10.0f-6, # Activation volume for coble diffusion creep
                p = 3f0,       # Grain size exponent
                n = 1f0,       # Stress exponent
                r = 0.7f0,     # Water fugacity exponent
                alf = 25f0,    # Melt factor
                ϕ_c = 1.0f-5,
                # x_ϕ_c = 5f0
            )
        ),

        disl = (
            dry = (
                A = 1.1f5,   # Preexponential
                Q = 530.0f3, # Activation energy
                V = 15.0f-6, # Activation volume
                p = 0f0,       # Grain size exponent
                n = 3.5f0,       # Stress exponent
                r = 0f0,       # Water fugacity exponent
                alf = 30f0,    # Melt factor
                ϕ_c = 1.0f-5,
                # x_ϕ_c = 1f0
            ),
            wet = (
                A = 1600.0f0,  # Preexponential
                Q = 520.0f3, # Activation energy
                V = 22.0f-6, # Activation volume                
                p = 0f0,       # Grain size exponent
                n = 3.5f0,     # Stress exponent
                r = 1.2f0,     # Water fugacity exponent
                alf = 30f0,    # Melt factor
                ϕ_c = 1.0f-5,
                # x_ϕ_c = 1f0
            )
        ),

        gbs = (
            lt1250 = (
                A = 6500.0f0,  # Preexponential for GBS-disl creep
                Q = 400.0f3, # Activation energy for GBS-disl creep
                V = 15.0f-6, # Activation volume
                p = 2f0,       # Grain size exponent
                n = 3.5f0,     # Stress exponent
                r = 0f0,       # Water fugacity exponent
                alf = 35f0,    # Melt factor
                ϕ_c = 1.0f-5,
                # x_ϕ_c = 25.0f-1
            ),
            gt1250 = (
                A = 4.7f10,  # Preexponential for GBS-disl creep
                Q = 600.0f3, # Activation energy for GBS-disl creep
                V = 15.0f-6, # Activation volume
                p = 2f0,       # Grain size exponent
                n = 3.5f0,     # Stress exponent
                r = 0f0,       # Water fugacity exponent
                alf = 35f0,    # Melt factor
                ϕ_c = 1.0f-5,
                # x_ϕ_c = 25.0f-1
            )
        )
    ),

    ch2o_o = 50,
    p_dep_calc = true,
    melt_enhancement = false
)


params_xfit_premelt = (
    α = 30f0, 
    T_η = 94f-2,
    γ = 5f0,
    B = 1f0,
    Tr = 1473f0,
    Pr = 1.5f0,
    η_r = 6.22f21,
    H = 462.5f3,
    V = 7.913f-6,
    M = 3f0,
    dg_r = 4f3
)

const default_params_HZK2011 = deepcopy(params_HZK2011)
const default_params_HK2003 = deepcopy(params_HK2003)
const default_params_xfit_premelt = deepcopy(params_xfit_premelt)
