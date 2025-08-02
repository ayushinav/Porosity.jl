function Base.show(
        io::IO, m::model) where {model <: MTModel{<:AbstractVector, <:AbstractVector}}
    println("1D ", typeof(m).name.name, " : ")
    println("Layer \t log(ρ)  h")
    println("____________________________")
    for i in eachindex(m.h)
        # println(resp_names_definition[k], " : ", getfield(m, k))
        m_ = round(m.m[i]; digits=3)
        h_ = round(m.h[i]; digits=3)
        println("$i \t $m_ \t $h_")
    end
    m_ = round(m.m[end]; digits=3)
    println("$(length(m.m)) \t $m_ \t ∞")
end
