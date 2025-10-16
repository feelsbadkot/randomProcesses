module ReturnTimeStatistics

using Peaks, Roots, StatsBase, Interpolations, Plots, LaTeXStrings, Measures

export plot_T_histogram, filling_factor, frequency_of_reversals

function plot_T_histogram(signal, time, δ=0.01)
    default(fontfamily="Computer Modern", titlefontsize=14, labelfontsize=12, legendfontsize=12, tickfontsize=10)
    indices, _ = findmaxima(signal)
    if length(indices) < 2
        return 0.0
    end
    T = diff(time[indices])
    Tmin, Tmax = extrema(T)
    return histogram(T, bins=range(Tmin, Tmax, step=δ), normalize=:probability, xlabel=L"T", ylabel=L"P(T)", legend=false, color="black", bottom_margin=3mm, left_margin=3mm)
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

end 