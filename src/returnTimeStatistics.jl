module ReturnTimeStatistics

using Peaks, Roots, BaseStats, Interpolations

function filling_factor(signal, time)
    indices, heights = findmaxima(signal)
    diff(time[indices])
    hist = fit(Histogram, data, nbins=bins)
    edges = hist.edges[1]
    bin_width = step(edges)
    Nt = count(!iszero, hist.weights)
    return Nt * bin_width
end 

function frequency_of_reversals(signal, time)
    f = linear_interpolation(time, X, extrapolation_bc=Line())
    find_zeros(f, time[0], time[end])
end

end 