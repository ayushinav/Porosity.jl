const model_names_definition = (
    T="Temperature (K)", ρ="Density (kg/m³)", P="Pressure (GPa)",
    ϕ="Porosity", dg="grain size(μm)", σ="Shear stress (GPa)",
    f="Frequency (Hz)", Ch2o_m="Water concentration in melt (ppm)",
    Ch2o_ol="Water concentration in olivine (ppm)",
    Ch2o_opx="Water concentration in orthopyroxene (ppm)",
    Ch2o_cpx="Water concentration in clinopyroxene (ppm)",
    Cco2_m="CO₂ concentration in melt (ppm)", T_solidus="Solidus Temperature (K)")

function Base.show(io::IO, m::model) where {model <: AbstractRockphyModel}
    println("Model : ", typeof(m).name.name)
    for k in propertynames(m)
        println(model_names_definition[k], " : ", getfield(m, k))
    end
end

const resp_names_definition = (
    T="Temperature (K)", σ="log₁₀ conductivity (S/m)", G="Elastic shear modulus (Pa)",
    K="Elastic bulk modulus (Pa)", Vp="Elastic P-wave velocity (m/s)",
    Vs="Elastic S-wave velocity (m/s)", ϵ_rate="Strain rate",
    η="Viscosity (Pa s)", J1="Real part of dynamic compliance (Pa⁻¹)",
    J2="Imaginary part of dynamic compliance (Pa⁻¹)", Qinv="Attenuation", M="Modulus (Pa)",
    V="Anelastic S-wave velocity : (m/s)", Vave="Frequency averaged S-wave velocity (m/s)")

function Base.show(io::IO, m::model) where {model <: AbstractRockphyResponse}
    println("Rock physics Response : ", typeof(m).name.name)
    for k in propertynames(m)
        println(resp_names_definition[k], " : ", getfield(m, k))
    end
end

function Base.show(io::IO, m::model) where {model <: two_phase_model}
    println("Two phase composition using ", m.mix, "\n")

    println("* ϕ₁ → ", 1 - first(m.ϕ), " : ", m.m1)
    println()
    println("* ϕ₂ → ", first(m.ϕ), " : ", m.m2)
end

function Base.show(io::IO, m::model) where {model <: multi_phase_model}
    println("N-phase composition using ", m.mix, "\n")

    fnames = propertynames(m)[2:(end - 1)]
    fnames = filter(f -> !isnothing(getfield(m, f)), fnames)

    for i in eachindex(fnames)
        sub = Char(0x2080 + i)
        println("* ϕ$sub → ", m.ϕ[i], " : ", getfield(m, fnames[i]))
        println()
    end
end

function Base.show(io::IO, m::model) where {model <: multi_rp_model}
    println("Multi rock physics model composed of \n")

    labels = (; cond="Conductivity model", elastic="Elastic model",
        visc="Viscosity model", anelastic="Anelastic model")

    for i in propertynames(m)
        if !isnothing(getfield(m, i))
            print("* ", getfield(labels, i), " : ")
            println(getfield(m, i))
        end
    end
end

function Base.show(io::IO, m::model) where {model <: multi_rp_response}
    println("Multi rock physics response composed of :  \n")

    labels = (; cond="Conductivity response", elastic="Elastic response",
        visc="Viscosity response", anelastic="Anelastic response")

    for i in propertynames(m)
        if !isnothing(getfield(m, i))
            print("* ", getfield(labels, i), " : ")
            println(getfield(m, i))
        end
    end
end
