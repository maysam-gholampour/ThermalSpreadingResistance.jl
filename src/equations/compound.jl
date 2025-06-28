function _Φ(::Compound, ξ, δₛ, hᶜ, 𝑘ₛ, 𝑘ₚ, δₚ)
    ϱ = (ξ + (hᶜ / 𝑘ₚ)) / (ξ - (hᶜ / 𝑘ₚ))
    κ = 𝑘ₚ / 𝑘ₛ
    α = (1.0 - κ) / (1.0 + κ)
    numerator_term_1 = α * exp(4ξ * δₛ) - exp(2ξ * δₛ)
    numerator_term_2 = ϱ * (exp(2ξ * (2δₛ + δₚ))-α*exp(2ξ * (δₛ + δₚ)))
    denominator_term_1 = α * exp(4ξ * δₛ) + exp(2ξ * δₛ)

    denominator_term_2 = ϱ * (exp(2ξ * (2δₛ + δₚ))+α*exp(2ξ * (δₛ + δₚ)))

    return (numerator_term_1 + numerator_term_2) / (denominator_term_1 + denominator_term_2)
end

function _get_Φ(::Compound, δₛ, hᶜ, 𝑘ₛ, 𝑘ₚ, δₚ)
    return ξ -> _Φ(Compound(), ξ, δₛ, hᶜ, 𝑘ₛ, 𝑘ₚ, δₚ)
end

function _Θ(::Compound, x, y, z, A₀, B₀, Aₘ, Aₙ, Aₘₙ, Bₘ, Bₙ, Bₘₙ, λ, δ, β)
    number_of_term = length(λ)
    sum_m = 0.0
    sum_n = 0.0
    sum_mn = 0.0
    @inbounds @fastmath @simd for i in 1:number_of_term
        sum_m += cos(λ[i] * x) * (Aₘ[i] * cosh(λ[i] * z) + Bₘ[i] * sinh(λ[i] * z))
        sum_n += cos(δ[i] * y) * (Aₙ[i] * cosh(δ[i] * z) + Bₙ[i] * sinh(δ[i] * z))
    end
    @inbounds @fastmath @simd for i in 1:number_of_term
        @inbounds @fastmath @simd for j in 1:number_of_term
            sum_mn += cos(λ[i] * x) *
                      cos(δ[j] * y) *
                      (Aₘₙ[i, j] * cosh(β[i, j] * z) + Bₘₙ[i, j] * sinh(β[i, j] * z))
        end
    end
    sss = A₀ + B₀ * z + sum_m + sum_n + sum_mn

    return sss
end

function _get_Θ(::Compound, A₀, B₀, Aₘ, Aₙ, Aₘₙ, Bₘ, Bₙ, Bₘₙ, λ, δ, β)
    return (x, y, z) -> _Θ(Isotropic(), x, y, z, A₀, B₀, Aₘ, Aₙ, Aₘₙ, Bₘ, Bₙ, Bₘₙ, λ, δ, β)
end

struct CompoundResults{T3 <: Function, T1 <: AbstractFloat}
    Θ::T3
    Θ_avg::T1
    Rₜ::T1
    R₁D::T1
    Rₛ::T1
end

function solve(::Compound, a, b, c, d, Q, 𝑘ₛ, δₛ, 𝑘ₚ, δₚ, hᶜ, Xᶜ, Yᶜ, number_of_term)
    A_b = a * b

    λ = zeros(number_of_term)
    δ = zeros(number_of_term)
    Aₘ = zeros(number_of_term)
    Aₙ = zeros(number_of_term)
    Bₘ = zeros(number_of_term)
    Bₙ = zeros(number_of_term)
    Aₘₙ = zeros((number_of_term, number_of_term))
    Bₘₙ = zeros((number_of_term, number_of_term))
    β = zeros((number_of_term, number_of_term))

    A₀ = (Q / (a * b)) * ((δₛ / 𝑘ₛ) + (1.0 / hᶜ))
    B₀ = -Q / (𝑘ₛ * a * b)

    Θ₁D = (Q / (a * b)) * ((δₚ / 𝑘ₚ) + (δₛ / 𝑘ₛ) + (1.0 / hᶜ))

    Φ = _get_Φ(Compound(), δₛ, hᶜ, 𝑘ₛ, 𝑘ₚ, δₚ)
    _calc_coefficients!(Aₘ, Aₙ, Aₘₙ, Bₘ, Bₙ, Bₘₙ, λ, δ, β, a, b, c, d, Q, 𝑘ₛ, Xᶜ, Yᶜ, Φ)

    Θ = _get_Θ(Compound(), A₀, B₀, Aₘ, Aₙ, Aₘₙ, Bₘ, Bₙ, Bₘₙ, λ, δ, β)  # Check Temperature distribution function

    Θ_avg = _Θ_avg(Aₘ, Aₙ, Aₘₙ, λ, δ, c, d, Xᶜ, Yᶜ, Θ₁D)
    Rₜ = Θ_avg / Q
    R₁D = (δₛ / (A_b * 𝑘ₛ)) + (δₚ / (A_b * 𝑘ₚ)) + (1.0 / (hᶜ * A_b))
    Rₛ = Rₜ - R₁D

    return CompoundResults(Θ, Θ_avg, Rₜ, R₁D, Rₛ)
end
