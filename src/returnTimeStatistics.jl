module ReturnTimeStatistics

using Peaks, Roots, StatsBase, Interpolations, Plots, LaTeXStrings, Measures

export plot_T_histogram, filling_factor, frequency_of_reversals, plot_filling_factor_map, plot_frequency_map

function plot_T_histogram(signal, time, δ=0.01)
    default(fontfamily="Computer Modern", titlefontsize=14, labelfontsize=12, legendfontsize=12, tickfontsize=10)
    indices, _ = findmaxima(signal)
    if length(indices) < 2
        return 0.0
    end
    T = diff(time[indices])
    Tmin, Tmax = extrema(T)
    return histogram(T, bins=range(Tmin, Tmax, step=δ), normalize=:probability, xlabel=L"T", ylabel=L"P(T)", 
                     legend=false, color="black", bottom_margin=3mm, left_margin=3mm)
end

function filling_factor(signal, time, δ=0.01)
    indices, _ = findmaxima(signal)
    if length(indices) < 2
        return 0.0
    end
    T = diff(time[indices])
    Tmin, Tmax = extrema(T)
    hist = fit(Histogram, T, range(Tmin, Tmax, step=δ))
    Nt = count(!iszero, hist.weights)
    return Nt * δ
end 

function frequency_of_reversals(signal, time)
    return length(find_zeros(linear_interpolation(time, signal, extrapolation_bc=Line()), time[1], time[end])) / findmax(time)[1]
end

function plot_filling_factor_map(σ_array, r_array, ff_array)
    default(fontfamily="Computer Modern", titlefontsize=14, labelfontsize=12, legendfontsize=12, tickfontsize=10)

    σ_grid = range(minimum(σ_array), maximum(σ_array), length=length(unique(σ_array)))
    r_grid = range(minimum(r_array), maximum(r_array), length=length(unique(r_array)))
    ff_grid = zeros(Float64, length(σ_grid), length(r_grid))

    k = 0
    for (i, σ) in enumerate(σ_grid)
        for (j, r) in enumerate(r_grid)
            k += 1
            if (σ == σ_array[k]) && (r == r_array[k])
                #println(true)
                ff_grid[i, j] = ff_array[k]
            end
        end
    end

    return contourf(σ_grid, r_grid, transpose(ff_grid),
                    levels=31, cmap=:hot, xlabel=L"\sigma", ylabel=L"r", colorbar=true, dpi=200)
end

end 