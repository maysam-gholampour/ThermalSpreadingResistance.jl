
@testset "Isotropic" begin
    a = 0.300
    b = 0.300
    c = 0.025
    d = 0.025
    Q = 10
    𝑘ₛ = 10
    δₛ = 10/1000
    hᶜ = 10
    Xᶜ = 0.090
    Yᶜ = 0.090
    number_terms=50

    Q2 = 15
    Xᶜ2 = 0.21
    Yᶜ2 = 0.21

    sol = solve(Isotropic(), a, b, c, d, Q, 𝑘ₛ, δₛ, hᶜ, Xᶜ, Yᶜ, number_terms)

    sol2 = solve(Isotropic(), a, b, c, d, Q2, 𝑘ₛ, δₛ, hᶜ, Xᶜ2, Yᶜ2, number_terms)

    @test sol.Θ_avg ≈ 47.7094 atol=1e-4
    @test sol.R₁D ≈ 1.1222 atol=1e-4
    @test sol.Rₛ ≈ 3.6487 atol=1e-4
    @test sol.Rₜ ≈ 4.7709 atol=1e-4

    @test sol.Θ(Xᶜ, Yᶜ, 0)+sol2.Θ(Xᶜ, Yᶜ, 0)+25 ≈ 84.9386 atol=1e-4
    @test sol.Θ(Xᶜ2, Yᶜ2, 0)+sol2.Θ(Xᶜ2, Yᶜ2, 0)+25 ≈ 108.38415 atol=1e-4
end
