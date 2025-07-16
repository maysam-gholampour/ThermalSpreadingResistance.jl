function _Φ(::Plate,ξ,δₛ,hᶜ,𝑘ₛ)
    
end
function _Φ(::Isotropic,ξ,δₛ,hᶜ,𝑘ₛ,𝑘ₚ=0.0,δₚ=0.0)
    return (ξ * sinh(ξ * δₛ) + (hᶜ / 𝑘ₛ) * cosh(ξ * δₛ)) / 
        (ξ * cosh(ξ * δₛ) + (hᶜ / 𝑘ₛ) * sinh(ξ * δₛ)) 
    
end
function _Φ(::Compound,ξ,δₛ,hᶜ,𝑘ₛ,𝑘ₚ,δₚ)
     ϱ = (ξ + (hᶜ / 𝑘ₚ)) / (ξ - (hᶜ / 𝑘ₚ))
    κ = 𝑘ₚ / 𝑘ₛ
    α = (1.0 - κ) / (1.0 + κ)
    numerator_term_1 = α * exp(4ξ * δₛ) - exp(2ξ * δₛ)
    numerator_term_2 = ϱ * (exp(2ξ * (2δₛ + δₚ))-α*exp(2ξ * (δₛ + δₚ)))
    denominator_term_1 = α * exp(4ξ * δₛ) + exp(2ξ * δₛ)

    denominator_term_2 = ϱ * (exp(2ξ * (2δₛ + δₚ))+α*exp(2ξ * (δₛ + δₚ)))

    return (numerator_term_1 + numerator_term_2) / (denominator_term_1 + denominator_term_2)
    
end
#a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ,Θ₁D,number_of_term
function _Θ_avg(plate,a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ,Θ₁D,number_of_term,𝑘ₚ::Float64=0.0, δₚ::Float64=0.0)
    
    sum_m = 0.0
    sum_n = 0.0
    sum_mn = 0.0
    @inbounds @fastmath @simd for i in 1:number_of_term
        λ = i * π / a
        Aₘ = 2 * Q * (
            sin(0.5 * (2Xᶜ + c) * λ) - sin(0.5 * (2Xᶜ - c) * λ)
            ) / (a * b * c * 𝑘ₛ * (λ ^ 2) * _Φ(plate,λ,δₛ,hᶜ,𝑘ₛ,𝑘ₚ, δₚ))
        
        δ = i * π / b
        Aₙ = 2 * Q * (
            sin(0.5 * (2Yᶜ + c) * δ) - sin(0.5 * (2Yᶜ - c) * δ)
            ) / (a * b * d * 𝑘ₛ * (δ ^ 2) * _Φ(plate,δ,δₛ,hᶜ,𝑘ₛ, 𝑘ₚ, δₚ))

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
                        (a * b * c * d * 𝑘ₛ * β * λ * δ * _Φ(plate,β,δₛ,hᶜ,𝑘ₛ,𝑘ₚ,δₚ))
            





                        
            sum_mn += Aₘₙ * cos(λ * Xᶜ) * cos(δ * Yᶜ) *
                        sin(0.5 * λ * c) * sin(0.5 * δ * d) / 
                        (λ * c * δ * d)
        end
    end
    return Θ₁D + 2sum_m + 2sum_n + 4sum_mn
end
# Ejemplo de cómo crear una instancia (llamar) a la función _Θ_avg
# Debes definir los valores de los parámetros antes de llamar a la función

# a = 1.0
# b = 1.0
# c = 0.5
# d = 0.5
# Q = 10.0
# 𝑘ₛ = 200.0
# Xᶜ = 0.25
# Yᶜ = 0.25
# Θ₁D = 0.0
# number_of_term = 10
# δₛ = 0.1
# hᶜ = 100.0

# # Llama a la función (esto devuelve el resultado)
# resultado = _Θ_avg(Isotropic(),a, b, c, d, Q, 𝑘ₛ,δₛ,hᶜ, Xᶜ, Yᶜ, Θ₁D, number_of_term)
# println("Resultado: ", resultado)
# # Ejemplo de cómo crear una instancia (llamar) a la función _Φ para Isotropic
# ξ = 1.0
# δₛ = 0.1
# hᶜ = 100.0
# 𝑘ₛ = 200.0


# resultado_Φ = _Φ(Isotropic(), ξ, δₛ, hᶜ, 𝑘ₛ)
# println("Resultado _Φ Isotropic: ", resultado_Φ)
