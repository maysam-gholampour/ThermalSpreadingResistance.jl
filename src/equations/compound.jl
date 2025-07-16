

function _calc_coefficients_vectorial!(i,a,b,c,d,Q,𝑘ₛ, δₛ,𝑘ₚ, δₚ, hᶜ,Xᶜ,Yᶜ)
    λ = i * π / a
    Aₘ = 2 * Q * (
            sin(0.5 * (2Xᶜ + c) * λ) - sin(0.5 * (2Xᶜ - c) * λ)
            ) / (a * b * c * 𝑘ₛ * (λ ^ 2) * _Φ(Compound(),λ,δₛ,hᶜ,𝑘ₛ,𝑘ₚ,δₚ))
    Bₘ = -_Φ(Compound(),λ,δₛ,hᶜ,𝑘ₛ,𝑘ₚ,δₚ) * Aₘ
    δ = i * π / b
    Aₙ = 2 * Q * (
        sin(0.5 * (2Yᶜ + c) * δ) - sin(0.5 * (2Yᶜ - c) * δ)
        ) / (a * b * d * 𝑘ₛ * (δ ^ 2) * _Φ(Compound(),δ,δₛ,hᶜ,𝑘ₛ,𝑘ₚ,δₚ))
    Bₙ = -_Φ(Compound(),δ,δₛ,hᶜ,𝑘ₛ,𝑘ₚ,δₚ) * Aₙ
    return Aₘ, Bₘ, Aₙ, Bₙ, λ, δ
end
function _calc_coefficients_matrix!(i,j,a,b,c,d,Q,𝑘ₛ,δₛ,𝑘ₚ, δₚ, hᶜ,Xᶜ,Yᶜ)
    λ = i * π / a
    δ = j * π / b
    β = √(λ^2 + δ^2)
    Aₘₙ = 16 * Q * cos(λ * Xᶜ) * sin(0.5 * λ * c) *
                cos(δ * Yᶜ) * sin(0.5 * δ * d) /
                (a * b * c * d * 𝑘ₛ * β * λ * δ * _Φ(Compound(),β,δₛ,hᶜ,𝑘ₛ,𝑘ₚ,δₚ))
    Bₘₙ = -_Φ(Compound(),β,δₛ,hᶜ,𝑘ₛ,𝑘ₚ,δₚ) * Aₘₙ
    return Aₘₙ, Bₘₙ, λ, δ, β
end


function _Θ(::Compound, x, y, z,A₀,B₀,a,b,c,d,Q,𝑘ₛ,δₛ,𝑘ₚ, δₚ, hᶜ,Xᶜ,Yᶜ,Θ₁D,number_of_term)
    
    sum_m = 0.0
    sum_n = 0.0
    sum_mn = 0.0
    @inbounds @fastmath @simd for i in 1:number_of_term

        Aₘ, Bₘ, Aₙ, Bₙ, λ, δ=_calc_coefficients_vectorial!(i,a,b,c,d,Q,𝑘ₛ,δₛ,𝑘ₚ, δₚ, hᶜ,Xᶜ,Yᶜ)


        sum_m += cos(λ * x) * (Aₘ * cosh(λ * z) + Bₘ * sinh(λ * z))
        sum_n += cos(δ * y) * (Aₙ * cosh(δ * z) + Bₙ * sinh(δ * z))
    end
    @inbounds @fastmath @simd for i in 1:number_of_term
        @inbounds @fastmath @simd for j in 1:number_of_term
        Aₘₙ, Bₘₙ, λ, δ, β=_calc_coefficients_matrix!(i,j,a,b,c,d,Q,𝑘ₛ,δₛ,𝑘ₚ, δₚ, hᶜ,Xᶜ,Yᶜ)

            sum_mn += cos(λ * x) *
                      cos(δ * y) *
                      (Aₘₙ * cosh(β * z) + Bₘₙ * sinh(β * z))
        end
    end
    sss = A₀ + B₀ * z + sum_m + sum_n + sum_mn

    return sss
end
#Θ = _get_Θ(Compound(), A₀,B₀,a,b,c,d,Q,𝑘ₛ,Xᶜ,Yᶜ,Φ,Θ₁D,number_of_term)
function _get_Θ(::Compound, A₀,B₀,a,b,c,d,Q,𝑘ₛ,δₛ,𝑘ₚ, δₚ, hᶜ,Xᶜ,Yᶜ,Θ₁D,number_of_term)
    return (x, y, z) -> _Θ(Compound(), x, y, z, A₀,B₀,a,b,c,d,Q,𝑘ₛ,δₛ,𝑘ₚ, δₚ, hᶜ,Xᶜ,Yᶜ,Θ₁D,number_of_term)
end

struct CompoundResults{T3 <: Function, T1 <: Real}
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

    
    #_calc_coefficients!(Aₘ, Aₙ, Aₘₙ, Bₘ, Bₙ, Bₘₙ, λ, δ, β, a, b, c, d, Q, 𝑘ₛ, Xᶜ, Yᶜ, Φ)

    Θ = _get_Θ(Compound(), A₀,B₀,a,b,c,d,Q,𝑘ₛ,δₛ,𝑘ₚ, δₚ, hᶜ,Xᶜ,Yᶜ,Θ₁D,number_of_term)  # Check Temperature distribution function
#plate,a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ,Θ₁D,number_of_term
    Θ_avg = _Θ_avg(Compound(),a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ,Θ₁D,number_of_term,𝑘ₚ, δₚ)
    Rₜ = Θ_avg / Q
    R₁D = (δₛ / (A_b * 𝑘ₛ)) + (δₚ / (A_b * 𝑘ₚ)) + (1.0 / (hᶜ * A_b))
    Rₛ = Rₜ - R₁D

    return CompoundResults(Θ, Θ_avg, Rₜ, R₁D, Rₛ)
end
