# MC_switch. First module for Monte-Carlo switching using Molecular Dynamics
module MC_switch

import Potentials
import Plots

function __init__()
  # Physical parameters
  global M = 1 # Mass of particles
  global kT = 1 # Boltzmann constant * Temperature

  # Numerical parameters
  global δt = 0.1 # Timestep length
  global notimesteps = 10000 # Number of timesteps to run for
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

# Gives the actual x-position of a displacement "dis" in a given well
function real_pos(well, dis)
  return dis + Potentials.minima[well]
end

# Gives the displacement in a given well given the real x-position
function well_dis(well, x_pos)
  return x_pos - Potentials.minima[well]
end

function simulate()
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

  for m = 1:notimesteps
    dis = well_dis(cur_well, x) # Displacement from the current well
    if m % switch_regularity == 0
      switch_attempts += 1

      ΔU = Potentials.U_shifted(x) - Potentials.U_shifted(real_pos(oth_well, dis)) # Energy change in making the switch
      if rand() <= min(1, exp(-ΔU))
        temp = cur_well
        cur_well = oth_well
        oth_well = temp
        x = real_pos(cur_well, dis)
      end

      if cur_well == 1
        noleft += 1
      end

      P_left = noleft / switch_attempts
      P_right = (switch_attempts - noleft) / switch_attempts
      push!(ΔF, -kT * log(P_left / P_right) + Potentials.U_shift)
    end

  x = x - (δt * Potentials.DU(x)) / M + sqrt((0.5 * kT * δt) / M) * (R + R_next)
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

end



end
