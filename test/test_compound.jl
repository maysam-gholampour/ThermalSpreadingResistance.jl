using Test
@testset "Compound" begin
    a = 0.05
    b = 0.05
    c = 0.025
    d = 0.025
    Q = 800.0
    𝑘ₛ = 200.0
    δₛ = 0.001
    𝑘ₚ = 5.0
    δₚ = 0.0001
    hᶜ = 1500.0

    Xᶜ = 0.5 * a
    Yᶜ = 0.5 * b

    sol = solve(Compound(), a, b, c, d, Q, 𝑘ₛ, δₛ, 𝑘ₚ, δₚ, hᶜ, Xᶜ, Yᶜ, 100)

    @test sol.Θ_avg ≈ 366.65039 atol=1e-4
    @test sol.R₁D ≈ 0.276666 atol=1e-5
    @test sol.Rₛ ≈ 0.181646 atol=1e-6
    @test sol.Rₜ ≈ 0.4583129 atol=1e-6
    @test sol.Θ(0, 0, 0) ≈ 103.89258364876 atol=1e-10
end
