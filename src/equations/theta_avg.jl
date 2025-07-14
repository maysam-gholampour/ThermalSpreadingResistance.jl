function _ΦΦ(::Isotropic,ξ,δₛ,hᶜ,𝑘ₛ)
    return (ξ * sinh(ξ * δₛ) + (hᶜ / 𝑘ₛ) * cosh(ξ * δₛ)) / 
        (ξ * cosh(ξ * δₛ) + (hᶜ / 𝑘ₛ) * sinh(ξ * δₛ)) 
    
end
function _Θ_avg(a,b,c,d,Q,𝑘ₛ,Xᶜ,Yᶜ,Θ₁D,number_of_term)
    
    sum_m = 0.0
    sum_n = 0.0
    sum_mn = 0.0
    @inbounds @fastmath @simd for i in 1:number_of_term
        λ = i * π / a
        Aₘ = 2 * Q * (
            sin(0.5 * (2Xᶜ + c) * λ) - sin(0.5 * (2Xᶜ - c) * λ)
            ) / (a * b * c * 𝑘ₛ * (λ ^ 2) * _Φ(λ,δₛ,hᶜ,𝑘ₛ))
        
        δ = i * π / b
        Aₙ = 2 * Q * (
            sin(0.5 * (2Yᶜ + c) * δ) - sin(0.5 * (2Yᶜ - c) * δ)
            ) / (a * b * d * 𝑘ₛ * (δ ^ 2) * _Φ(δ,δₛ,hᶜ,𝑘ₛ))

        sum_m += Aₘ * cos(λ * Xᶜ) * sin(0.5 * λ * c) / (λ * c)
        sum_n += Aₙ * cos(δ * Yᶜ) * sin(0.5 * δ * d) / (δ * d)
    end
    @inbounds @fastmath @simd for i in 1:number_of_term
        λ = i * π / a
        @inbounds @fastmath @simd for j in 1:number_of_term
            δ = j * π / b
            β = √(λ^2 + δ^2)
            Aₘₙ = 16 * Q * cos(λ * Xᶜ) * sin(0.5 * λ * c) * 
                        cos(δ * Yᶜ) * sin(0.5 * δ * d) / 
                        (a * b * c * d * 𝑘ₛ * β * λ * δ * _Φ(β,δₛ,hᶜ,𝑘ₛ))
            





                        
            sum_mn += Aₘₙ * cos(λ * Xᶜ) * cos(δ * Yᶜ) *
                        sin(0.5 * λ * c) * sin(0.5 * δ * d) / 
                        (λ * c * δ * d)
        end
    end
    return Θ₁D + 2sum_m + 2sum_n + 4sum_mn
end