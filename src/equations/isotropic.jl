



function _calc_coefficients_vectorial!(i,a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ)
    λ = i * π / a
    Aₘ = 2 * Q * (
            sin(0.5 * (2Xᶜ + c) * λ) - sin(0.5 * (2Xᶜ - c) * λ)
            ) / (a * b * c * 𝑘ₛ * (λ ^ 2) * _Φ(Isotropic(),λ,δₛ,hᶜ,𝑘ₛ))
    Bₘ = -_Φ(Isotropic(),λ,δₛ,hᶜ,𝑘ₛ) * Aₘ
    δ = i * π / b
    Aₙ = 2 * Q * (
        sin(0.5 * (2Yᶜ + c) * δ) - sin(0.5 * (2Yᶜ - c) * δ)
        ) / (a * b * d * 𝑘ₛ * (δ ^ 2) * _Φ(Isotropic(),δ,δₛ,hᶜ,𝑘ₛ))
    Bₙ = -_Φ(Isotropic(),δ,δₛ,hᶜ,𝑘ₛ) * Aₙ
    return Aₘ, Bₘ, Aₙ, Bₙ, λ, δ
end

function _calc_coefficients_matrix!(i,j,a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ)
    λ = i * π / a
    δ = j * π / b
    β = √(λ^2 + δ^2)
    Aₘₙ = 16 * Q * cos(λ * Xᶜ) * sin(0.5 * λ * c) * 
                cos(δ * Yᶜ) * sin(0.5 * δ * d) / 
                (a * b * c * d * 𝑘ₛ * β * λ * δ * _Φ(Isotropic(),β,δₛ,hᶜ,𝑘ₛ))
    Bₘₙ = -_Φ(Isotropic(),β,δₛ,hᶜ,𝑘ₛ) * Aₘₙ
    return Aₘₙ, Bₘₙ, λ, δ, β
end
function _Θ(::Isotropic,x,y,z,A₀,B₀,a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ,Θ₁D,number_of_term) #a,b,c,d,Q,𝑘ₛ,Xᶜ,Yᶜ,Φ,Θ₁D,number_of_term
    
    sum_m = 0.0
    sum_n = 0.0
    sum_mn = 0.0
    @inbounds  @fastmath @simd for i in 1:number_of_term

        Aₘ, Bₘ, Aₙ, Bₙ, λ, δ=_calc_coefficients_vectorial!(i,a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ)
        
        
        sum_m += cos(λ * x) * (Aₘ * cosh(λ * z) + Bₘ * sinh(λ * z))
        sum_n += cos(δ * y) * (Aₙ * cosh(δ * z) + Bₙ * sinh(δ * z))
    end
    @inbounds @fastmath @simd  for i in 1:number_of_term
        
        @inbounds @fastmath @simd  for j in 1:number_of_term
        Aₘₙ, Bₘₙ, λ, δ, β=_calc_coefficients_matrix!(i,j,a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ)    
        sum_mn += cos(λ * x) * cos(δ * y) * 
                (Aₘₙ * cosh(β * z) + Bₘₙ * sinh(β * z))
        end
    end
     sss= A₀ + B₀ * z + sum_m + sum_n + sum_mn
    return sss
end

function _get_Θ(::Isotropic,A₀,B₀,a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ,Θ₁D,number_of_term)
    return (x,y,z) -> _Θ(Isotropic(),x,y,z,A₀,B₀,a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ,Θ₁D,number_of_term)
end

struct IsoResults{T1<:Function,T2<:Real} 
    Θ::T1
    Θ_avg::T2
    Rₜ::T2
    R₁D::T2
    Rₛ::T2
end

function solve(::Isotropic,a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ,number_of_term)
    A_b = a * b

    A₀ = (Q / (a * b)) * ((δₛ / 𝑘ₛ) + (1.0 / hᶜ))
    B₀ = -Q / (𝑘ₛ * a * b)
    Θ₁D = A₀

    #Φ = _get_Φ(Isotropic(),δₛ,hᶜ,𝑘ₛ)
    Θ = _get_Θ(Isotropic(),A₀,B₀,a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ,Θ₁D,number_of_term)
    Θ_avg = _Θ_avg(Isotropic(),a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ,Θ₁D,number_of_term)
    Rₜ = Θ_avg / Q
    R₁D = (δₛ / (A_b * 𝑘ₛ)) + (1.0 / (hᶜ * A_b))
    Rₛ = Rₜ - R₁D
    return IsoResults(Θ,Θ_avg,Rₜ,R₁D,Rₛ)
end

    