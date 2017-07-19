# MC_switch. First module for Monte-Carlo switching using Molecular Dynamics
module MC_switch

import Potentials
import Plots

function __init__()
  # Physical parameters
  global M = 1 # Mass of particles
  global kT = 1 # Boltzmann constant * Temperature

  # Numerical parameters
  global δt = 0.01 # Timestep length
  global notimesteps = 1000000 # Number of timesteps to run for
  global switch_regularity = 10 # Number of timesteps between each switch attempt
  global x_min = -4
  global x_max = 4.4
  global start_well = 1 # Which well to start in, 1 is left, 2 is right
end

function plot()
  x = linspace(x_min, x_max, 1000)
  y = [Potentials.U_shifted(p) for p in x]
  Plots.plot(x, y)
end

function simulate(potential_name)
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
  switch_attempts = 0

  cur_well = start_well # Current well, 1 for left, 2 for right
  if cur_well == 1
    oth_well = 2
  else
    oth_well = 1
  end


  R = randn()
  R_next = randn()

  # x = 2.9272649686148524
  # println(DU(x))
  # x = x - (δt * DU(x)) / M + sqrt((0.5 * kT * δt) / M) * (R + R_next)
  # println(x)


  for m = 1:notimesteps
    if isnan(x)
      return("You're fuggered now")
    end
    # println(string("x = ", string(x)))
    # println(string("cur_well = ", string(cur_well)))

    if m % switch_regularity == 0
      switch_attempts += 1

      # Recalibrate the wells, if a particle has managed to cross over the barrier
      if (cur_well == 1 && x > 0) || (cur_well == 2 && x < 0)
        cur_well, oth_well = oth_well, cur_well
      end
      
      dis = well_dis(cur_well, x) # Displacement from the current well
      ΔU = U_shifted(real_pos(oth_well, dis)) - U_shifted(x) # Energy change in making the switch
      # println(" ")
      # println(ΔU)
      # println(exp(-ΔU))
      # println(" ")
      if rand() <= min(1, exp(-ΔU))
        cur_well, oth_well = oth_well, cur_well
        x = real_pos(cur_well, dis)
      end

      if cur_well == 1
        noleft += 1
      end

      P_left = noleft / switch_attempts
      P_right = (switch_attempts - noleft) / switch_attempts
      # println([P_left, P_right])
      push!(ΔF, -kT * log(P_left / P_right) + U_shift)
    end

    x = x - (δt * DU(x)) / M + sqrt((0.5 * kT * δt) / M) * (R + R_next)

    R = R_next
    R_next = randn()
  end
  # Clean up ΔF
  for j = [1:length(ΔF);]
    if ΔF[j] == Inf
      ΔF[j] = 0
    end
  end
  plt = Plots.plot([1:switch_attempts;], ΔF)
  Plots.display(plt)
  println(last(ΔF))
end

function exact_sol(potential_name)
  potential_arr = Potentials.potential_selector(potential_name)
  U = potential_arr[3]
  integrand(x) = exp(- U(x) / kT)
  P_left = quadgk(integrand, -1000, 0)
  P_right = quadgk(integrand, 0, 1000)
  println([P_left, P_right])
  println(-kT * log(P_left[1] / P_right[1]))
end

end
