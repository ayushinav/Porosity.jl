"""
`normal_dist(resp::AbstractVector, err_resp::AbstractVector)`: returns a multivariate normal distribution centered around `resp` with std dev `err_resp`
"""
function normal_dist(resp::AbstractVector, err_resp::AbstractVector)
    MultivariateNormal(resp, Diagonal(err_resp))
end

function normal_dist(resp::AbstractVector, err_resp::AbstractMatrix)
    MultivariateNormal(resp, err_resp)
end

"""
`uniform_dist(resp::AbstractVector, err_resp::AbstractVector)`: returns a uniform normal distribution ∈ [`resp`-`err_resp`/2, `resp`+`err_resp`/2]
"""
function uniform_dist(resp::AbstractVector, err_resp::AbstractVector)
    Product([Uniform(resp[k] - err_resp[k] / 2, resp[k] + err_resp[k] / 2)
             for k in eachindex(resp)])
end
