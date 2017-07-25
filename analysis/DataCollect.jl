# DataCollect.jl
#= Module for conducting large numbers of MC switch simulations, calculating free energy differences
 at varying temperatures and timestep lengths =#

module DataCollect

import MC_switch
import Auto_Correlate

# Runs simulations and stores result in a CSV
function simulate(potential_name, data_dir, timerange, timeinc, temprange, tempinc, maxtimesteps; filetag="")
    min_time = timerange[1]
    max_time = timerange[2]
    max_temp = temprange[2]


    start_time = min_time # Only relevant when starting from an existiing file MAKE THIS MORE ELEGANT
    CATCHUP = false

    filename = string(potential_name, "_", filetag, "_data.csv")
    fullfilename = joinpath(data_dir, filename)
    if isfile(fullfilename)
        endline = readdlm(fullfilename, '\t')[end, :]
        start_time = endline[2] + timeinc
        if start_time > timerange[2]
            min_temp = endline[1] + tempinc
            start_time = timerange[1]
            if min_temp > temprange[2]
                return "All simulations completed"
            end 
        else
            min_temp = endline[1]
            CATCHUP = true
        end

    else
        min_temp = temprange[1]
    end
    
    interval_length = max(1, round(Int64, maxtimesteps / 10000))

    outputfilename = joinpath(data_dir, string(potential_name, "_output.csv"))
    for kT = [min_temp:tempinc:max_temp;]
        exact_val = MC_switch.exact_sol(potential_name, kT)
        for δt = [min_time:timeinc:max_time;]
            if (CATCHUP) && (δt < start_time)
                continue
            end

            output_data = []
            for i = [1:3;]
                try
                    MC_switch.simulate(potential_name, maxtimesteps, δt, kT, outputfilename)
                catch e
                    println("Caught error $e")
                    continue
                end
                push!(output_data, Auto_Correlate.calculate(outputfilename, 0.5, 0.01, interval_length=interval_length))
            end

            result_file = open(fullfilename, "a+")

            if length(output_data) == 0
                println("No finite values reached")
                writedlm(result_file, [kT, δt, NaN, NaN, exact_val]')
            else
                println(length(output_data))
                println(output_data)
                meanval = mean(output_data[j][1] for j = 1:length(output_data))
                stderr = sqrt(sum(output_data[j][2]^2 for j = 1:length(output_data))) / (length(output_data))
                writedlm(result_file, [kT, δt, meanval, stderr, exact_val]')
            end

            close(result_file)

            if (CATCHUP) && (δt == max_time)
                CATCHUP = false
            end
        end
    end
end


end
