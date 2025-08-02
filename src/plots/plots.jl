## response plots

"""
    plot_response!(axs, vars, resp::response; errs=zero(resp), plt_type=:plot,
        kwargs...)

plots on `axs`, a vector of `axis` corresponding parameters of response `resp`
parameterized on `vars`

## Arguments

  - `axs` : vector of `axis` to plot on
  - `vars` : parameter for forward sol, eg. ω for MT, ab/2 for DC resistivity
  - `resp` : response to be plotted

## Keyword Arguments

  - `errs` : errors in response, should be of the same type as `resp`, defaults to `zero(resp)`

  - `plt_type` : type of plot, available options are :

      + `:scatter` : data points are plotted as points
      + `:plot` : data points are connected in a curve
      + `:errors` : plot errors along with the data points
        Defaults to `:plot`
  - `kwargs` : keyword arguments to be splatted into the `Makie` plotting routine

## Usage

Checkout relevant documentation
"""
function plot_response!(axs, vars, resp::response; errs=zero(resp), plt_type=:plot,
        kwargs...) where {response <: AbstractGeophyResponse}
    k = fieldnames(response)

    if plt_type === :plot
        for i in eachindex(axs)
            lines!(axs[i], vars, getproperty(resp, k[i]); kwargs...)
            xscale, yscale = get_scales(response, Val{k[i]}())

            axs[i].xscale = xscale
            axs[i].yscale = yscale

            xlabel, ylabel = get_labels(response, Val{k[i]}())

            axs[i].xlabel = xlabel
            axs[i].ylabel = ylabel
        end

    elseif plt_type === :scatter
        for i in eachindex(axs)
            scatter!(axs[i], vars, getproperty(resp, k[i]); kwargs...)
            xscale, yscale = get_scales(response, Val{k[i]}())

            axs[i].xscale = xscale
            axs[i].yscale = yscale

            xlabel, ylabel = get_labels(response, Val{k[i]}())

            axs[i].xlabel = xlabel
            axs[i].ylabel = ylabel
        end
    elseif plt_type === :errors
        for i in eachindex(axs)
            errorbars!(axs[i], vars, getproperty(resp, k[i]),
                getproperty(errs, k[i]) ./ 2; kwargs...)
            xscale, yscale = get_scales(response, Val{k[i]}())

            axs[i].xscale = xscale
            axs[i].yscale = yscale

            xlabel, ylabel = get_labels(response, Val{k[i]}())

            axs[i].xlabel = xlabel
            axs[i].ylabel = ylabel
        end
    end
end

"""
    plot_response(vars, resp::response; errs=zero(resp), plt_type=:plot,
        kwargs...)

returns a figure and `axis` with plots of corresponding parameters of response `resp`
parameterized on `vars`

## Arguments

  - `vars` : parameter for forward sol, eg. ω for MT, ab/2 for DC resistivity
  - `resp` : response to be plotted

## Keyword Arguments

  - `errs` : errors in response, should be of the same type as `resp`, defaults to `zero(resp)`

  - `plt_type` : type of plot, available options are :

      + `:scatter` : data points are plotted as points
      + `:plot` : data points are connected in a curve
      + `:errors` : plot errors along with the data points
        Defaults to `:plot`
  - `kwargs` : keyword arguments to be splatted into the `Makie` plotting routine

## Usage

Checkout relevant documentation
"""
function plot_response(vars, resp::response; errs=zero(resp), plt_type=:plot,
        kwargs...) where {response <: AbstractGeophyResponse}
    k = fieldnames(response)
    f = Figure()
    axs = [Axis(f[i, 1]) for i in eachindex(k)]

    plot_response!(axs, vars, resp; errs=errs, plt_type=plt_type, kwargs...)

    return f, axs
end

## model plots
"""
    plot_model!(ax, model::m_type; half_space_thickness=1.25 * sum(model.h),
        kwargs...)

plots on `ax` the geophysical model parameterized by `model`

## Arguments

  - `ax` : `axis` on which models are plotted
  - `model` : geophysical model to be plotted

## Keyword Arguments

  - `half_space_thickness` : thickness of last layer (extent to which the last layer
    aka the half space should go), defaults to 1.25 × total model thickness

  - `kwargs` : keyword arguments to be splatted into the `Makie` plotting routine

## Usage

Checkout relevant documentation
"""
function plot_model!(ax, model::m_type; half_space_thickness=1.25 * sum(model.h),
        kwargs...) where {m_type}
    # ax = f.content[1]
    m = model.m
    h = model.h

    m_vec = 10 .^ [m[1], m...]
    h_v = cumsum(h)
    h_vec = [1.0f-2, h_v..., half_space_thickness]

    stairs!(ax, m_vec, h_vec; step=:post, kwargs...)
    ax.yreversed = true
    nothing
end

"""
    plot_model(model::m_type; half_space_thickness=1.25 * sum(model.h),
        kwargs...)

returns a `figure` and an `axis` with the plot of geophysical model parameterized by `model`

## Arguments

  - `ax` : `axis` on which models are plotted
  - `model` : geophysical model to be plotted

## Keyword Arguments

  - `half_space_thickness` : thickness of last layer (extent to which the last layer
    aka the half space should go), defaults to 1.25 × total model thickness

  - `kwargs` : keyword arguments to be splatted into the `Makie` plotting routine

## Usage

Checkout relevant documentation
"""
function plot_model(
        model::m_type; half_space_thickness=1.25 * sum(model.h), kwargs...) where {m_type}
    fig = Figure()
    ax = Axis(fig[1, 1])

    plot_model!(ax, model; half_space_thickness=half_space_thickness, kwargs...)
    ax.yreversed = true

    fig, ax
end
