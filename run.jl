include("src/lorenzSystem.jl")
include("src/returnTimeStatistics.jl")

import .LorenzSystem as LS
import .ReturnTimeStatistics as RTS

using HDF5

struct LorenzData{T<:Real}
    σ::T
    r::T
    b::T
    solution::Dict
    filling_factor::T
    frequency_of_reversals::T
end

function parse_range(str::String)
    """
    Парсит строку вида "start:step:stop" или "start,stop,length"
    """
    if ':' in str
        # Формат start:step:stop
        parts = split(str, ':')
        start = parse(Float64, parts[1])
        step = parse(Float64, parts[2])
        stop = parse(Float64, parts[3])
        return range(start, step=step, stop=stop)
    elseif ',' in str
        # Формат start,stop,length
        parts = split(str, ',')
        start = parse(Float64, parts[1])
        stop = parse(Float64, parts[2])
        length = parse(Int, parts[3])
        return range(start, stop=stop, length=length)
    else
        error("Неверный формат диапазона. Используйте 'start:step:stop' или 'start,stop,length'")
    end
end

function main(σ_range_str::String, r_range_str::String, 
              b::Float64=8/3, 
              tspan_start::Float64=0.0, tspan_stop::Float64=2500.0,
              u0::Vector{Float64}=[0.0, 1.0, 0.0],
              output_filename::String="lorenz_results.h5")
    
    # Парсим диапазоны из строк
    σ_range = parse_range(σ_range_str)
    r_range = parse_range(r_range_str)
    tspan = (tspan_start, tspan_stop)
    
    # Создаем HDF5 файл
    #h5open(output_filename, "cw") do file
    open(output_filename, "w") do file
        
        # Предварительное создание массивов для результатов
        n_total = length(σ_range) * length(r_range)
        sigmas = Vector{Float64}(undef, n_total)
        rs = Vector{Float64}(undef, n_total)
        filling_factors = Vector{Float64}(undef, n_total)
        reversal_frequencies = Vector{Float64}(undef, n_total)
        
        idx = 1
        for σ in σ_range
            for r in r_range
                p = (σ, r, b)
                
                # Решаем систему
                solution = LS.solve_lorenz(u0, p, tspan)
                solution_dict = LS.solution_to_dict(solution)
                
                # Вычисляем статистики
                filling = RTS.filling_factor(solution_dict["X"], solution_dict["t"])
                frequency = RTS.frequency_of_reversals(solution_dict["X"], solution_dict["t"])
                
                # Сохраняем в массивы
                
                println(file, "$idx,$σ,$r,$filling,$frequency")
                println("$idx \t $σ \t $r \t $filling \t $frequency")

                idx += 1
            end
        end
        
        # Записываем результаты в файл
        #file["sigma"] = sigmas
        #file["r"] = rs
        #file["filling_factor"] = filling_factors
        #file["reversal_frequency"] = reversal_frequencies
        #file["b"] = [b]  # Сохраняем как массив из одного элемента
        #
        ## Добавляем метаданные
        #attrs(file)["sigma_range"] = σ_range_str
        #attrs(file)["r_range"] = r_range_str
        #attrs(file)["time_span"] = "$tspan"
        #attrs(file)["initial_conditions"] = string(u0)
        #attrs(file)["output_file"] = output_filename
        
        println("Results are saved in $output_filename")
    end
end

# Запуск
σ_range_str = ARGS[1]
r_range_str = ARGS[2]

b = 8 / 3

tspan_start = length(ARGS) >= 3 ? parse(Float64, ARGS[3]) : 0.0
tspan_stop = length(ARGS) >= 4 ? parse(Float64, ARGS[4]) : 3000.0

u0 = [0.0, 1.0, 0.0]

output_file = length(ARGS) >= 5 ? ARGS[5] : "lorenz_results.h5"

main(σ_range_str, r_range_str, b, tspan_start, tspan_stop, [0.0, 1.0, 0.0], output_file)