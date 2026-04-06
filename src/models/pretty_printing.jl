# COV_EXCL_START
const model_names_definition = (
    T="Temperature (K)", ρ="Density (kg/m³)", P="Pressure (GPa)",
    cond="Conductivity (S/m)", ϕ="Porosity", dg="grain size(μm)", σ="Shear stress (GPa)",
    f="Frequency (Hz)", Ch2o_m="Water concentration in melt (ppm)",
    Ch2o_ol="Water concentration in olivine (ppm)",
    Ch2o_opx="Water concentration in orthopyroxene (ppm)",
    Ch2o_cpx="Water concentration in clinopyroxene (ppm)",
    Cco2_m="CO₂ concentration in melt (ppm)", T_solidus="Solidus Temperature (K)")

function Base.show(
        io::IO, ::MIME"text/plain", m::model) where {model <: AbstractRockphyModel}
    println(io, "Model : ", typeof(m).name.name)
    for k in propertynames(m)
        println(io, model_names_definition[k], " : ", getfield(m, k))
    end
end

const resp_names_definition = (
    T="Temperature (K)", σ="log₁₀ conductivity (S/m)", G="Elastic shear modulus (Pa)",
    K="Elastic bulk modulus (Pa)", Vp="Elastic P-wave velocity (m/s)",
    Vs="Elastic S-wave velocity (m/s)", ϵ_rate="Strain rate",
    η="Viscosity (Pa s)", J1="Real part of dynamic compliance (Pa⁻¹)",
    J2="Imaginary part of dynamic compliance (Pa⁻¹)", Qinv="Attenuation", M="Modulus (Pa)",
    V="Anelastic S-wave velocity : (m/s)", Vave="Frequency averaged S-wave velocity (m/s)")

function Base.show(
        io::IO, ::MIME"text/plain", m::model) where {model <: AbstractRockphyResponse}
    println(io, "Rock physics Response : ", typeof(m).name.name)
    for k in propertynames(m)
        println(io, resp_names_definition[k], " : ", getfield(m, k))
    end
end

function Base.show(io::IO, ::MIME"text/plain", m::model) where {model <: two_phase_model}
    println(io, "Two phase composition using ", m.mix, "\n")
    ϕ_ = @. 1 - m.ϕ

    println(io, "*    m₁ (solid phase) : ")
    show(io, MIME"text/plain"(), m.m1)
    println(io, "ϕ : ", ϕ_)
    println(io)
    println(io, "*    m₂ (liquid phase) : ")
    show(io, MIME"text/plain"(), m.m1)
    println(io, "ϕ : ", m.ϕ)
end

function Base.show(io::IO, ::MIME"text/plain", m::model) where {model <: multi_phase_model}
    println(io, "multi-phase composition using ", m.mix)

    fnames = propertynames(m)[2:(end - 1)]
    fnames = filter(f -> !isnothing(getfield(m, f)), fnames)

    for i in eachindex(fnames)
        sub = Char(0x2080 + i)
        println(io)
        println(io, "*    m$sub : ")
        show(io, MIME"text/plain"(), getfield(m, fnames[i]))
        println(io, "ϕ : ", m.ϕ[i])
    end
end

function Base.show(io::IO, ::MIME"text/plain", m::model) where {model <: multi_rp_model}
    println(io, "Multi rock physics model composed of \n")

    labels = (; cond="Conductivity model", elastic="Elastic model",
        visc="Viscosity model", anelastic="Anelastic model")

    for i in propertynames(m)
        if !isnothing(getfield(m, i))
            print("* ", getfield(labels, i), " : ")
            println(getfield(m, i))
        end
    end
end

function Base.show(io::IO, ::MIME"text/plain", m::model) where {model <: tune_rp_modelType}
    println(io, "Tuning rock physics with function list :")
    println(io, m.fn_list)
    println(io, "to obtain the rock physics model of type ", m.model)
end

function Base.show(io::IO, ::MIME"text/plain", m::model) where {model <: multi_rp_response}
    println(io, "Multi rock physics response composed of :  \n")

    labels = (; cond="Conductivity response", elastic="Elastic response",
        visc="Viscosity response", anelastic="Anelastic response")

    for i in propertynames(m)
        if !isnothing(getfield(m, i))
            print("* ", getfield(labels, i), " : ")
            println(getfield(m, i))
        end
    end
end
# COV_EXCL_STOP
