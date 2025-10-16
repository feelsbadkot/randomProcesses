module LorenzSystem 

using DifferentialEquations, Plots, LaTeXStrings, Measures

export solve_lorenz, solution_to_dict, plot_time_series, plot_phase_portrait

const PLOT_INDEX_MAP = Dict(
    "X" => [1], "1" => [1], 1 => [1],
    "Y" => [2], "2" => [2], 2 => [2], 
    "Z" => [3], "3" => [3], 3 => [3],
    "3D" => [4], "4" => [4], 4 => [4], "XYZ3D" => [4],
    
    "XY" => [1,2], "12" => [1,2],
    "XZ" => [1,3], "13" => [1,3],
    "YZ" => [2,3], "23" => [2,3],
    "X3D" => [1,4], "14" => [1,4],
    "Y3D" => [2,4], "24" => [2,4],
    "Z3D" => [3,4], "34" => [3,4],
    
    "XYZ" => [1,2,3], "123" => [1,2,3],
    "XY3D" => [1,2,4], "124" => [1,2,4],
    "XZ3D" => [1,3,4], "134" => [1,3,4],
    "YZ3D" => [2,3,4], "234" => [2,3,4],
    
    "ALL" => [1,2,3,4], "all" => [1,2,3,4], "1234" => [1,2,3,4]
)

function lorenz!(du, u, p, t)
    σ, r, b = p
    du[1] = σ * (u[2] - u[1])
    du[2] = u[1] * (r - u[3]) - u[2]
    du[3] = u[1] * u[2] - b * u[3]
end

function solve_lorenz(u0, p, tspan)
    prob = ODEProblem(lorenz!, u0, tspan, p)
    t_output = range(tspan[1], tspan[2], step=0.01)
    return solve(prob, saveat=t_output) 
end

function solution_to_dict(sol)
    return Dict("t" => sol.t, "X" => sol[1, :], "Y" => sol[2, :], "Z" => sol[3, :])
end

function plot_time_series(sol, index, limits=extrema(sol.t))
    default(fontfamily="Computer Modern", titlefontsize=14, labelfontsize=12, legendfontsize=12, tickfontsize=10)
    time_series = (plot(sol, plotdensity = 10000, lw = 1.25, idxs = (0, 1), xlabel = L"t", ylabel = L"X(t)", xlims=limits, legend = false, color="black", bottom_margin=3mm, left_margin=3mm, right_margin=3mm),
                   plot(sol, plotdensity = 10000, lw = 1.25, idxs = (0, 2), xlabel = L"t", ylabel = L"Y(t)", xlims=limits, legend = false, color="black", bottom_margin=3mm, left_margin=3mm, right_margin=3mm), 
                   plot(sol, plotdensity = 10000, lw = 1.25, idxs = (0, 3), xlabel = L"t", ylabel = L"Z(t)", xlims=limits, legend = false, color="black", bottom_margin=3mm, left_margin=3mm, right_margin=3mm),
                   plot(sol, plotdensity = 10000, lw = 1.25, xlabel = L"t", xlims=limits, legend = false, bottom_margin=3mm, left_margin=3mm, right_margin=3mm))
    plot_indices = get(PLOT_INDEX_MAP, index, [1, 2, 3, 4]) 
    return time_series[plot_indices]
end

function plot_phase_portrait(sol, index)
    default(fontfamily="Computer Modern", titlefontsize=14, labelfontsize=12, legendfontsize=12, tickfontsize=10)
    
    phase_plots = (plot(sol, plotdensity=10000, lw=1.25, idxs=(1, 2), xlabel=L"X", ylabel=L"Y", legend=false, color="black", bottom_margin=3mm, left_margin=3mm, right_margin=3mm),
                   plot(sol, plotdensity=10000, lw=1.25, idxs=(1, 3), xlabel=L"X", ylabel=L"Z", legend=false, color="black", bottom_margin=3mm, left_margin=3mm, right_margin=3mm),
                   plot(sol, plotdensity=10000, lw=1.25, idxs=(2, 3), xlabel=L"Y", ylabel=L"Z", legend=false, color="black", bottom_margin=3mm, left_margin=3mm, right_margin=3mm),
                   plot(sol, plotdensity=10000, lw=1, idxs=(1, 2, 3), xlabel=L"X", ylabel=L"Y", zlabel=L"Z", legend=false, color="black"))
    
    plot_indices = get(PLOT_INDEX_MAP, index, [1, 2, 3, 4])  
    return phase_plots[plot_indices]
end

end 