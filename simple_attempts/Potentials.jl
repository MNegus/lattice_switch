module Potentials


function __init__()
  # Potential specific parameters
  global minima = [-2, 2.4] # First element is location of left minimum, other is right minumum
  global U_shift = 4.4 # Amount the right mimumum has been shifted up
end

# External potential function
function U(x)
  if x <= -1
    return 5 * (x + 2)^2
  elseif x <= 1.2
    return 10 - 5 * x^2
  else
    return 5 * (x - 2.4)^2 - 4.4
  end
end

# Potential function shifted
function U_shifted(x)
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
function DU(x)
  if x <= -1
    return 10 * (x + 2)
  elseif x <= 1.2
    return -10 * x
  else
    return 10 * (x - 2.4)
  end
end

end
