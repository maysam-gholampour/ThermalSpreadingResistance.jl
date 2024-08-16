function _Θ_avg(Aₘ,Aₙ,Aₘₙ,λ,δ,c,d,Xᶜ,Yᶜ,Θ₁D)
    number_of_term = length(λ)
    sum_m = 0.0
    sum_n = 0.0
    sum_mn = 0.0
    @inbounds @fastmath @simd for i in 1:number_of_term
        sum_m += Aₘ[i] * cos(λ[i] * Xᶜ) * sin(0.5 * λ[i] * c) / (λ[i] * c)
        sum_n += Aₙ[i] * cos(δ[i] * Yᶜ) * sin(0.5 * δ[i] * d) / (δ[i] * d)
    end
    @inbounds @fastmath @simd for i in 1:number_of_term
        @inbounds @fastmath @simd for j in 1:number_of_term
            sum_mn += Aₘₙ[i,j] * cos(λ[i] * Xᶜ) * cos(δ[j] * Yᶜ) *
                        sin(0.5 * λ[i] * c) * sin(0.5 * δ[j] * d) / 
                        (λ[i] * c * δ[j] * d)
        end
    end
    return Θ₁D + 2sum_m + 2sum_n + 4sum_mn
end