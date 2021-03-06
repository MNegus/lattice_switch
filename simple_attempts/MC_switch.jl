# MC_switch. First module for Monte-Carlo switching using Molecular Dynamics
module MC_switch

import Potentials
# import Plots
import Cuba
# Plots.pyplot()

function __init__()
  # Physical parameters
  global M = 1 # Mass of particles

  # Numerical parameters
  global x_min = -4
  global x_max = 1.5
  global start_well = 1 # Which well to start in, 1 is left, 2 is right
end

# function plot()
#   x = linspace(x_min, x_max, 1000)
#   y = [Potentials.U_shifted(p) for p in x]
#   Plots.plot(x, y)
# end

function simulate(potential_name, maxtimesteps, δt, kT; switch_regularity=10)
  potential_arr = Potentials.potential_selector(potential_name)
  minima = potential_arr[1]
  U_shift = potential_arr[2]
  U = potential_arr[3]
  U_shifted = potential_arr[4]
  DU = potential_arr[5]
  real_pos(well, dis) = dis + minima[well]
  well_dis(well, x_pos) = x_pos - minima[well]

  x = real_pos(start_well, 0) # Starting position of the walker


  ΔF = [] # Array to store difference in free energy after every switch attempt
  noleft = 0 # Number of switch attempts the walker has been in left well

  cur_well = start_well # Current well, 1 for left, 2 for right
  if cur_well == 1
    oth_well = 2
  else
    oth_well = 1
  end

  NO_SHIFTS = 0

  R = randn()
  R_next = randn()
  
  for timestep = 1:maxtimesteps
    if cur_well == 1
      noleft += 1
    end


    # BAOAB step
    x = x - (δt * DU(x)) / M + sqrt((0.5 * kT * δt) / M) * (R + R_next)
    R = R_next
    R_next = randn()
    if isnan(x)
      println(string("Infinite value of x reached at timestep ", timestep))   
      throw(DomainError())
    end

    # Recalibrate the wells, if a particle has managed to cross over the barrier
    if (cur_well == 1 && x > 0) || (cur_well == 2 && x < 0)
      cur_well, oth_well = oth_well, cur_well
    end


    
    if timestep % switch_regularity == 0
      dis = well_dis(cur_well, x) # Displacement from the current well
      ΔU = U_shifted(real_pos(oth_well, dis)) - U_shifted(x) # Energy change in making the switch
      # println(" ")
      # println(string("m = ", m))
      # println(string("cur_well = ", cur_well))
      # println(string("dis = ", dis))
      # println(string("x = ", x))
      # println(string("ΔU = ", ΔU))
      # println(string("exp(-ΔU) = ", exp(-ΔU)))
      # println(" ")
      if rand() <= min(1, exp(-ΔU))
        cur_well, oth_well = oth_well, cur_well
        x = real_pos(cur_well, dis)
      end

      P_left = noleft / timestep
      P_right = (timestep - noleft) / timestep
      push!(ΔF, -kT * log(P_left / P_right) + U_shift)
    end    
  end

  # Clean up ΔF
  for j = [1:length(ΔF);]
    if isinf(ΔF[j])
      ΔF[j] = 0
    end
  end

  # writedlm(filename, ΔF, "\t")
  # plt = Plots.plot([1:switch_attempts;], ΔF)
  # Plots.savefig(string(potential_name, "_ΔF.png"))
  return convert(Array{Float64, 1}, ΔF)
end

function exact_sol(potential_name, kT)
  potential_arr = Potentials.potential_selector(potential_name)
  int_func(x) = exp(- potential_arr[3](x) / kT)
  P_left = Cuba.vegas((t, f) -> f[1] = int_func(-t[1] / (1 - t[1])) / ((1 - t[1])^2), reltol=10.0^(-7))
  P_right = Cuba.vegas((t, f) -> f[1] = int_func(t[1] / (1 - t[1])) / ((1 - t[1])^2), reltol=10.0^(-7))
return -kT * log(P_left[1][1] / P_right[1][1])
end

end
