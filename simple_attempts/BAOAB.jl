# BAOAB.jl
# Module for implementing BAOAB steps, both using the BAOAB limit method and regular BAOAB

module BAOAB

# Calculates the next x values using BAOAB limit method
function limit(x, δt, M, ∇U, kT, R, R_next)
    return x - δt * ∇U(x) / M + sqrt(0.5 * kT * δt / M) * (R + R_next)
end

# Calculates the next x and p values using the regular BAOAB method
function standard(p, x, δt, M, ∇U, c1, c2, c3, R_next)
    p_half = p - 0.5 * δt * ∇U(x)
    x_half = x + 0.5 * δt * p_half / M
    p_half_bar = c1 * p_half + c3 * sqrt(M) * R_next
    return (p_half_bar - 0.5 * δt * ∇U(x), x_half + 0.5 * δt * p_half_bar / M)
end
end
