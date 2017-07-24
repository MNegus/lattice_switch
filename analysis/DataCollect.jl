# DataCollect.jl
#= Module for conducting large numbers of MC switch simulations, calculating free energy differences
 at varying temperatures and timestep lengths =#

module DataCollect

import MC_switch
import Auto_Correlate

# Runs simulations and stores result in a CSV
function simulate(potential_name, timerange, timeinc, temprange, tempinc, maxtimesteps)
    min_time = timerange[1]
    max_time = timerange[2]
    max_temp = temprange[2]


    start_time = min_time # Only relevant when starting from an existiing file MAKE THIS MORE ELEGANT

    data_dir = "/home/michael/Documents/lattice_switch/analysis/" # Directory data stored in
    filename = string(potential_name, "_data.csv")
    fullfilename = joinpath(data_dir, filename)
    if isfile(fullfilename)
        endline = readdlm(fullfilename, ',')[end, :]
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

    result_file = open(fullfilename, "a+")
    
    interval_length = max(1, round(Int64, maxtimesteps / 10000))

    outputfilename = joinpath(data_dir, string(potential_name, "_output.csv"))
    for kT = [min_temp:tempinc:max_temp;]
        for δt = [min_time:timeinc:max_time]
            println(kT, δt)
            if (CATCHUP) && (δt < start_time)
                continue
            end

            output_data = []
            for i = [1:3;]
                MC_switch.simulate(potential_name, maxtimesteps, δt, kT, outputfilename)
                println(Auto_Correlate.calculate(outputfilename, 0.5, 0.01, interval_length=interval_length))
            end
            print(output_data)
            if (CATCHUP) && (δt == max_time)
                CATCHUP = false
            end
        end
    end

    close(result_file)

end


end