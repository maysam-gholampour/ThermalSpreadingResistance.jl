function phi_zeta(zeta, t1, t2, k1, k2, h)
    kappa = k2 / k1
    alpha = (1 - kappa) / (1 + kappa)
    eg= (zeta+h/k2) / (zeta-h/k2)

    num = (alpha * exp(4*zeta * t1) - exp(2*zeta * t1)) 
    + eg*(exp(2 * zeta * (2*t1+t2)) - alpha*exp(2 * zeta * (t1+t2)))
    den = (alpha * exp(4*zeta * t1) + exp(2*zeta * t1)) 
    + eg *(exp(2 * zeta * (2*t1+t2)) + alpha*exp(2 * zeta * (t1+t2)))
    
    if den == 0
      return NaN # Handle division by zero
    else
      phi= num / den
      return phi
    end
end
  
function theta_1D_avg( k1, k2, h, t1, t2)
    temperature = h*(t1 / k1) + h*(t2 / k2) + 1
    return temperature
end

# example usage
#=zeta_val = 1.5
t1_val = 0.1
t2_val = 0.2
k1_val = 10.0
k2_val = 20.0
h_val = 100.0

phi_at_zeta = phi_zeta(zeta_val, t1_val, t2_val, k1_val, k2_val, h_val)
println("Value phi(zetaa)  = $zeta_val: $phi_at_zeta")


theta_avg = theta_1D_avg( k1_val, k2_val, h_val, t1_val, t2_val)
println("Value of θ̄₁  for compound system: $theta_avg")
=#