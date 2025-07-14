using Test

@testset "Isotropic" begin
    a = 0.05
    b = 0.05
    c = 0.025
    d = 0.025
    Q = 800.0
    𝑘ₛ = 200.0
    δₛ = 0.001
    hᶜ = 150000.0

    Xᶜ = 0.5 * a
    Yᶜ = 0.5 * b

    sol = solve(Isotropic(), a, b, c, d, Q, 𝑘ₛ, δₛ, hᶜ, Xᶜ, Yᶜ, 100)

    @test sol.Θ_avg ≈ 13.5885 atol=1e-4
    @test sol.R₁D ≈ 0.00466 atol=1e-5
    @test sol.Rₛ ≈ 0.012319 atol=1e-6
    @test sol.Rₜ ≈ 0.016985 atol=1e-6
    @test sol.Θ(0, 0, 0) ≈ 3.51175e-5 atol=1e-10
end





