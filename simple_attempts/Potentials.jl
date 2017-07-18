module Potentials


# Function which returns the relevant functions for a given named potential
function potential_selector(potential_name)
  if potential_name == "KT_NOTES"
    minima = [-2, 2.4]
    U_shift = 4.4
    return (minima, U_shift, KT_NOTES_U, KT_NOTES_U_shifted, KT_NOTES_DU)
  elseif potential_name == "QUARTIC"
    minima = [-1.220997215942, 1.5989977937]
    U_shift = 2.8256360858458973
    return (minima, U_shift, QUARTIC_U, QUARTIC_U_shifted, QUARTIC_DU)
  end
end


#=
 Potential function that appeared in the kinetic theory notes
=#
# External potential function
function KT_NOTES_U(x)
  if x <= -1
    return 5 * (x + 2)^2
  elseif x <= 1.2
    return 10 - 5 * x^2
  else
    return 5 * (x - 2.4)^2 - 4.4
  end
end

# Potential function shifted
function KT_NOTES_U_shifted(x)
  if x <= -1
    return 5 * (x + 2)^2
  elseif x <= 0
    return 10 - 5 * x^2
  elseif x <= 1.2
    return 14.4 - 5 * x^2
  else
    return 5 * (x - 2.4)^2
  end
end

# Derivative of potential function
function KT_NOTES_DU(x)
  if x <= -1
    return 10 * (x + 2)
  elseif x <= 1.2
    return -10 * x
  else
    return 10 * (x - 2.4)
  end
end



#=
 Quartic potential function
 =#
# Potential function
function QUARTIC_U(x)
  x_0 = -0.126000192586256
  return (x + x_0)^4 - 4 * (x + x_0)^2 - (x + x_0) + 2.618555980765
end

# Potential function shifted
function QUARTIC_U_shifted(x)
  x_0 = -0.126000192586256
  if x < 0
    return (x + x_0)^4 - 4 * (x + x_0)^2 - (x + x_0) + 2.618555980765
  else
    return (x + x_0)^4 - 4 * (x + x_0)^2 - (x + x_0) + 5.444192066610897
  end
end

# Derivative of potential function
function QUARTIC_DU(x)
  x_0 = -0.126000192586256
  return 4 * (x + x_0)^3 - 8 * (x + x_0) - 1
end

end
