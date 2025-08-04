using ThermalSpreadingResistance
using Test
using JET
using Enzyme

@testset "Stability with JET" begin
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
    number_of_terms = 100

    #Isotropic stability tests
    sol = @inferred ThermalSpreadingResistance.solve(
        ThermalSpreadingResistance.Isotropic(), a, b, c, d, Q, 𝑘ₛ, δₛ, hᶜ, Xᶜ, Yᶜ, number_of_terms)
    report_Θ_avg = @report_opt ThermalSpreadingResistance._Θ_avg(
        Isotropic(), a, b, c, d, Q, 𝑘ₛ, δₛ, hᶜ, Xᶜ, Yᶜ, 0.01, number_of_terms)
    report_get_Θ = @report_opt ThermalSpreadingResistance._get_Θ(
        ThermalSpreadingResistance.Isotropic(), 4.3, 3.2, a, b, c,
        d, Q, 𝑘ₛ, δₛ, hᶜ, Xᶜ, Yᶜ, 0.01, number_of_terms)
    report_solve = @report_opt sol.Rₛ

    num_issues_Θ_avg = length(JET.get_reports(report_Θ_avg))
    num_issues_get_Θ = length(JET.get_reports(report_get_Θ))
    num_issues_solve = length(JET.get_reports(report_solve))

    @test isa(sol, ThermalSpreadingResistance.IsoResults)
    @test num_issues_Θ_avg ≈ 0
    @test num_issues_get_Θ ≈ 0
    @test num_issues_solve ≈ 0

    #Compound stability tests
    𝑘ₚ = 5.0
    δₚ = 0.0001

    sol_compound = @inferred ThermalSpreadingResistance.solve(
        ThermalSpreadingResistance.Compound(), a, b, c, d,
        Q, 𝑘ₛ, δₛ, 𝑘ₚ, δₚ, hᶜ, Xᶜ, Yᶜ, number_of_terms)

    report_Θ_avg_compound = @report_opt ThermalSpreadingResistance._Θ_avg(
        Compound(), a, b, c, d, Q, 𝑘ₛ, δₛ, hᶜ, Xᶜ, Yᶜ, 0.01, number_of_terms, 𝑘ₚ, δₚ)
    report_get_Θ_compound = @report_opt ThermalSpreadingResistance._get_Θ(
        ThermalSpreadingResistance.Compound(), 4.3, 3.2, a, b, c,
        d, Q, 𝑘ₛ, δₛ, 𝑘ₚ, δₚ, hᶜ, Xᶜ, Yᶜ, 0.01, 100)

    report_solve_compound = @report_opt sol_compound.Rₛ

    num_issues_Θ_avg_compound = length(JET.get_reports(report_Θ_avg_compound))
    num_issues_get_Θ_compound = length(JET.get_reports(report_get_Θ_compound))
    num_issues_solve_compound = length(JET.get_reports(report_solve_compound))

    @test isa(sol_compound, ThermalSpreadingResistance.CompoundResults)
    @test num_issues_Θ_avg_compound ≈ 0
    @test num_issues_get_Θ_compound ≈ 0
    @test num_issues_solve_compound ≈ 0
end
