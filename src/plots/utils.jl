get_scales(::Type{<:MTResponse}, ::Val{:ρₐ}) = log10, log10
get_scales(::Type{<:MTResponse}, ::Val{:ϕ}) = log10, identity

get_labels(::Type{<:MTResponse}, ::Val{:ρₐ}) = "T(s)", "ρₐ (Ωm)"
get_labels(::Type{<:MTResponse}, ::Val{:ϕ}) = "T(s)", "ϕ (ᴼ)"

get_scales(::Type{<:MTModel}) = log10, log10
get_labels(::Type{<:MTModel}) = "ρ (Ωm)", "h(m)"
