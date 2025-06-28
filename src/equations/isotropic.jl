
function _get_Φ(::Isotropic, δₛ, hᶜ, 𝑘ₛ)
    _Φ(::Isotropic,
        ξ,
        δₛ,
        hᶜ,
        𝑘ₛ) = (ξ * sinh(ξ * δₛ) + (hᶜ / 𝑘ₛ) * cosh(ξ * δₛ)) /
              (ξ * cosh(ξ * δₛ) + (hᶜ / 𝑘ₛ) * sinh(ξ * δₛ))
    return ξ -> _Φ(Isotropic(), ξ, δₛ, hᶜ, 𝑘ₛ)
end

function _Θ(::Isotropic, x, y, z, A₀, B₀, Aₘ, Aₙ, Aₘₙ, Bₘ, Bₙ, Bₘₙ, λ, δ, β)
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
    return A₀ + B₀ * z + sum_m + sum_n + sum_mn
end

function _get_Θ(::Isotropic, A₀, B₀, Aₘ, Aₙ, Aₘₙ, Bₘ, Bₙ, Bₘₙ, λ, δ, β)
    return (x, y, z) -> _Θ(Isotropic(), x, y, z, A₀, B₀, Aₘ, Aₙ, Aₘₙ, Bₘ, Bₙ, Bₘₙ, λ, δ, β)
end

struct IsoResults{T1 <: Function, T2 <: AbstractFloat}
    Θ::T1
    Θ_avg::T2
    Rₜ::T2
    R₁D::T2
    Rₛ::T2
end

function solve(::Isotropic, a, b, c, d, Q, 𝑘ₛ, δₛ, hᶜ, Xᶜ, Yᶜ, number_of_term)
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
    Θ₁D = A₀

    Φ = _get_Φ(Isotropic(), δₛ, hᶜ, 𝑘ₛ)
    _calc_coefficients!(Aₘ, Aₙ, Aₘₙ, Bₘ, Bₙ, Bₘₙ, λ, δ, β, a, b, c, d, Q, 𝑘ₛ, Xᶜ, Yᶜ, Φ)

    Θ = _get_Θ(Isotropic(), A₀, B₀, Aₘ, Aₙ, Aₘₙ, Bₘ, Bₙ, Bₘₙ, λ, δ, β)
    Θ_avg = _Θ_avg(Aₘ, Aₙ, Aₘₙ, λ, δ, c, d, Xᶜ, Yᶜ, Θ₁D)
    Rₜ = Θ_avg / Q
    R₁D = (δₛ / (A_b * 𝑘ₛ)) + (1.0 / (hᶜ * A_b))
    Rₛ = Rₜ - R₁D
    return IsoResults(Θ, Θ_avg, Rₜ, R₁D, Rₛ)
end
