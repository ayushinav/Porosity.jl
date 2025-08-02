# Working with `mtpy`

`MT.jl` allows users to take advantage of the [mtpy](https://github.com/MTgeophysics/mtpy-v2) through the python calls offered by [PyCall.jl](https://github.com/JuliaPy/PyCall.jl). Installing `Pycall.jl` and making `mtpy` work can be slightly tricky. As `mtpy` gets updated, users would want to use the latest updated versions. To make things easier, we do not install `mtpy` via conda or python, and not get involved in the messy virtual environments, we simply clone `mtpy` from the github [repo](https://github.com/MTgeophysics/mtpy-v2).

*Add compat dependency for PyCall*

We provide a walktthrough to help get things started and get to the point to use this package. We begin by importing different packages.

```julia
using LinearAlgebra
using MT
using PyCall
using Statistics
using Plots
using JLD2
```

To import `mtpy`, we need to provide its path, relative paths work as well. In our case, the repo is cloned at `../../mtpy_v2/mtpy`.

We now have access to all the functions provided by `mtpy`. You can now use the package any way you want. The following lines of code demonstrate how we take a `.h5` file, which is also created by `mtpy` (check relevant tutorials).

```julia
np = pyimport("numpy")

filename = abspath(joinpath(dirname(@__FILE__), "../../../../mtpy_v2/mtpy"))
(path, name) = dirname(filename), basename(filename)

@pyimport imp
(file, filename, data) = imp.find_module(name, [path])
mtpy = imp.load_module(name, file, filename, data)
```

```julia
mc = mtpy.MTCollection()

mc.open_collection("../../../data_for_mtpy_demo/tf_collection2.h5")

df = mc.master_dataframe
lats = 1 .* mc.master_dataframe["latitude"]
lons = 1 .* mc.master_dataframe["longitude"]

idx_sort = 1 .+ np.argsort(lats);
n_stations = length(idx_sort);
st_keys0 = np.array(mc.master_dataframe["survey"])
st_keys = st_keys0[idx_sort] # station keys, we call st 0 as unknown_survey, st 1 as unknown_survey_002
```

At the time of writing this, I was not able to completely figure out how to export apparent resistivity and phase data along with the corresponding error bars. In the next section, we calculate these parameters using julia and demonstrate compatibility with `mtpy` functions and objects.

```julia
tf_arr = [];
z_arr = [];
zerr_arr = [];
f_arr = [];
for i in 1:n_stations
    tfi = mc.get_tf(; tf_id="P", survey=st_keys[i])
    # tf_arr.append(tfi);
    push!(tf_arr, tfi)
    push!(z_arr, np.array(tfi.impedance))
    push!(zerr_arr, np.array(tfi.impedance_error))
    # zerr_arr.append(np.array(tfi.impedance_error));
    push!(f_arr, np.array(tfi.frequency))
end

ρₐ_all = [abs.(z_arr[i]) .^ 2 ./ (2 * π .* f_arr[i]) for i in eachindex(z_arr)]
ϕ_all = [angle.(z_arr[i]) .* 180 / π for i in eachindex(z_arr)]

for i in eachindex(ϕ_all)
    ϕ_all[i][ϕ_all[i] .< 0.0] .= ϕ_all[i][ϕ_all[i] .< 0.0] .+ 180.0
end

ρₐerr = [2 .* zerr_arr[i] .* abs.(z_arr[i]) ./ (2 * π .* f_arr[i]) for i in 1:n_stations]
ϕerr = [];
for i in 1:n_stations
    r = zerr_arr[i] ./ abs2.(z_arr[i])
    err_r = zero(r)
    err_r[r .> 1] .= NaN
    err_r[r .<= 1] .= 180 / π .* asin.(zerr_arr[i][r .<= 1] ./ abs2.(z_arr[i][r .<= 1]))
    push!(ϕerr, err_r)
end

plts = [];
using Plots
plt3 = plot();
for i in 1:n_stations
    scatter!(plt3, [lons[i]], [lats[i]]; markersize=2, label=false)
end

anim = @animate for i in 1:n_stations
    plt1 = scatter(1 ./ f_arr[i], ρₐ_all[i][:, 1, 2]; yerr=(ρₐerr[i][:, 1, 2]),
        scale=:log10, label="XY", markersize=3)
    scatter!(plt1, 1 ./ f_arr[i], ρₐ_all[i][:, 2, 1];
        yerr=(ρₐerr[i][:, 2, 1]), scale=:log10, label="YX", markersize=3)
    plot!(plt1; ylim=(1, 1e4), size=(1200, 500), ylabel="ρₐ(Ωm)", xlabel="T(s)")

    plt2 = scatter(1 ./ f_arr[i], ϕ_all[i][:, 1, 2]; yerr=(ϕerr[i][:, 1, 2]),
        xscale=:log10, label="XY", markersize=3)
    scatter!(plt2, 1 ./ f_arr[i], ϕ_all[i][:, 2, 1]; yerr=(ϕerr[i][:, 2, 1]),
        xscale=:log10, label="YX", markersize=3)
    plot!(plt2; ylim=(0, 90), size=(1200, 500), ylabel="ϕ(°)", xlabel="T(s)")

    plt = plot(plt1, plt2; layout=(1, 2), margin=5Plots.mm,
        legend=:outertopright, size=(800, 350), title="$i")
end

gif(anim; fps=0.5)
```

Similar plots can also be made with `MT.jl`. We just need to create the corresponding `MTResponse`s and plot them.

```julia
anim = @animate for i in 1:n_stations
    xy_resp = MTResponse(ρₐ_all[i][:, 1, 2], ϕ_all[i][:, 1, 2])
    xy_err_resp = MTResponse(ρₐerr[i][:, 1, 2], ϕerr[i][:, 1, 2])

    yx_resp = MTResponse(ρₐ_all[i][:, 2, 1], ϕ_all[i][:, 2, 1])
    yx_err_resp = MTResponse(ρₐerr[i][:, 2, 1], ϕerr[i][:, 2, 1])

    plt = prepare_plot(xy_resp, 2π .* f_arr[i], xy_err_resp;
        label="XY", markershape=:circle, maerkersize=3)
    prepare_plot!(plt, yx_resp, 2π .* f_arr[i], yx_err_resp;
        label="YX", markershape=:circle, markersize=3)
    plot!(plt[1]; ylim=(1, 1e5))
    plot!(plt[2]; ylim=(0, 90))

    plot_response(plt; layout=(1, 2), margin=5Plots.mm,
        legend=:outertopright, size=(800, 350), title="$i")
end

gif(anim; fps=0.5)
```
