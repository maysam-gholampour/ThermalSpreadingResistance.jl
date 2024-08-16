function _Φ(::Compound,ξ,δₛ,hᶜ,𝑘ₛ,𝑘ₚ,δₚ)
    ϱ = (ξ + (hᶜ / 𝑘ₚ)) / (ξ - (hᶜ / 𝑘ₚ))
    κ = 𝑘ₚ / 𝑘ₛ
    α = (1.0 - κ) / (1.0 + κ)
    numerator_term_1 = α * exp(4ξ * δₛ) - exp(2ξ * δₛ)
    numerator_term_2 = ϱ * (1.0 - α) * exp(2ξ * (δₛ + δₚ))
    denominator_term_1 = α * exp(4ξ * δₛ) + exp(2ξ * δₛ)
    denominator_term_2 = ϱ * (1.0 + α) * exp(2ξ * (δₛ + δₚ))
    
    return (numerator_term_1 + numerator_term_2) / 
            (denominator_term_1 + denominator_term_2)
end

function _get_Φ(::Compound,δₛ,hᶜ,𝑘ₛ,𝑘ₚ,δₚ)
    return ξ -> _Φ(Compound(),ξ,δₛ,hᶜ,𝑘ₛ,𝑘ₚ,δₚ)
end

struct ComponudResults{T1<:AbstractFloat} 
    Θ_avg::T1
    Rₜ::T1
    R₁D::T1
    Rₛ::T1
end

function solve(::Compound,a,b,c,d,Q,𝑘ₛ,δₛ,𝑘ₚ,δₚ,hᶜ,number_of_term)
    Xᶜ = 0.5 * a
    Yᶜ = 0.5 * b
    A_b = a * b

    λ = zeros(number_of_term)
    δ = zeros(number_of_term)
    Aₘ = zeros(number_of_term)
    Aₙ = zeros(number_of_term)
    Bₘ = zeros(number_of_term)
    Bₙ = zeros(number_of_term)
    Aₘₙ = zeros((number_of_term,number_of_term))
    Bₘₙ = zeros((number_of_term,number_of_term))
    β = zeros((number_of_term,number_of_term))


    Θ₁D = (Q / (a * b)) * ((δₚ / 𝑘ₚ) + (δₛ / 𝑘ₛ) + (1.0 / hᶜ))

    Φ = _get_Φ(Compound(),δₛ,hᶜ,𝑘ₛ,𝑘ₚ,δₚ)
    _calc_coefficients!(Aₘ,Aₙ,Aₘₙ,Bₘ,Bₙ,Bₘₙ,λ,δ,β,a,b,c,d,Q,𝑘ₛ,Xᶜ,Yᶜ,Φ)

    Θ_avg = _Θ_avg(Aₘ,Aₙ,Aₘₙ,λ,δ,c,d,Xᶜ,Yᶜ,Θ₁D)
    Rₜ = Θ_avg / Q
    R₁D = (δₛ / (A_b * 𝑘ₛ)) + (δₚ / (A_b * 𝑘ₚ)) + (1.0 / (hᶜ * A_b))
    Rₛ = Rₜ - R₁D
    return ComponudResults(Θ_avg,Rₜ,R₁D,Rₛ)
end

