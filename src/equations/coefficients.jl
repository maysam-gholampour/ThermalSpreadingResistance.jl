function _calc_coefficients!(Aₘ,Aₙ,Aₘₙ,Bₘ,Bₙ,Bₘₙ,λ,δ,β,a,b,c,d,Q,𝑘ₛ,Xᶜ,Yᶜ,Φ)
    number_of_term = length(λ)
    @inbounds @fastmath @simd for i in 1:number_of_term
        λ[i] = i * π / a
        Aₘ[i] = 2 * Q * (
            sin(0.5 * (2Xᶜ + c) * λ[i]) - sin(0.5 * (2Xᶜ - c) * λ[i])
            ) / (a * b * c * 𝑘ₛ * (λ[i] ^ 2) * Φ(λ[i]))
        Bₘ[i] = -Φ(λ[i]) * Aₘ[i]
        δ[i] = i * π / b
        Aₙ[i] = 2 * Q * (
            sin(0.5 * (2Yᶜ + c) * δ[i]) - sin(0.5 * (2Yᶜ - c) * δ[i])
            ) / (a * b * d * 𝑘ₛ * (δ[i] ^ 2) * Φ(δ[i]))
        Bₙ[i] = -Φ(δ[i]) * Aₙ[i]
    end
    @inbounds @fastmath @simd for i in 1:number_of_term
        @inbounds @fastmath @simd for j in 1:number_of_term
            β[i,j] = √(λ[i]^2 + δ[j]^2)
            Aₘₙ[i,j] = 16 * Q * cos(λ[i] * Xᶜ) * sin(0.5 * λ[i] * c) * 
                        cos(δ[j] * Yᶜ) * sin(0.5 * δ[j] * d) / 
                        (a * b * c * d * 𝑘ₛ * β[i,j] * λ[i] * δ[j] * Φ(β[i,j]))
            Bₘₙ[i,j] = -Φ(β[i,j]) * Aₘₙ[i,j]
        end
    end
    return nothing
end
