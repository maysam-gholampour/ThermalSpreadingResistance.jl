
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
    hᶜ = 150000.0

    Xᶜ = 0.5 * a
    Yᶜ = 0.5 * b
    
    sol = solve(Compound(),a,b,c,d,Q,𝑘ₛ,δₛ,𝑘ₚ,δₚ,hᶜ,Xᶜ,Yᶜ,100)

    @test sol.Θ_avg ≈ 42.4596 atol=1e-4
    @test sol.R₁D ≈ 0.01266 atol=1e-5
    @test sol.Rₛ ≈ 0.040407 atol=1e-6
    @test sol.Rₜ ≈ 0.053074 atol=1e-6
end



