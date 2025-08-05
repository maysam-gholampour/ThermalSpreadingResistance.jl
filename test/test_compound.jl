using Test
@testset "Compound" begin
    a = 0.300
    b = 0.300
    c = 0.025
    d = 0.025
    Q = 10.0
    𝑘ₛ = 350.0
    δₛ = 2/1000
    𝑘ₚ = 10.0
    δₚ = 10/1000
    hᶜ = 10.0
    Xᶜ = 0.090
    Yᶜ = 0.090
    number_terms=50

    Q2 = 15
    Xᶜ2 = 0.21
    Yᶜ2 = 0.21

    sol = solve(Compound(), a, b, c, d, Q, 𝑘ₛ, δₛ, 𝑘ₚ, δₚ, hᶜ, Xᶜ, Yᶜ, number_terms)

    sol2 = solve(Compound(), a, b, c, d, Q2, 𝑘ₛ, δₛ, 𝑘ₚ, δₚ, hᶜ, Xᶜ2, Yᶜ2, number_terms)

    @test sol.Θ_avg ≈ 16.18995 atol=1e-4
    @test sol.R₁D ≈ 1.12228 atol=1e-5
    @test sol.Rₛ ≈ 0.4967 atol=1e-4
    @test sol.Rₜ ≈ 1.61899 atol=1e-5

    @test sol.Θ(Xᶜ, Yᶜ, 0)+sol2.Θ(Xᶜ, Yᶜ, 0)+25 ≈ 56.2702 atol=1e-4
    @test sol.Θ(Xᶜ2, Yᶜ2, 0)+sol2.Θ(Xᶜ2, Yᶜ2, 0)+25 ≈ 59.66968 atol=1e-4
end
