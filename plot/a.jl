include("../src/lorenzSystem.jl")
include("../src/returnTimeStatistics.jl")

import .LorenzSystem as LS
import .ReturnTimeStatistics as RTS

using CSV, DataFrames, Plots

data = CSV.read("results/output_6.csv", DataFrame, header=false)
sigmas = data[!, "Column2"]
rs = data[!, "Column3"]
fillings = data[!, "Column4"];

p = RTS.plot_filling_factor_map(sigmas, rs, fillings)
png(p, "a")