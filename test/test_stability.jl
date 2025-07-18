using ThermalSpreadingResistance
using Test
using JET 
using Enzyme

@testset "Stability with JET" begin
        

    x = [
    0.008,  # a
    0.008,  # b
    0.0025, # c
    0.0025, # d
    20.0,   # Q
    300,    # 𝑘ₛ
    0.003,  # δₛ
    50.0,   # hᶜ
    0.004,  # Xᶜ
    0.004   # Yᶜ
    ]
    number_of_term = 100

        #Isotropic stability tests
        report_Θ_avg= @report_opt ThermalSpreadingResistance._Θ_avg(Isotropic(), x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], 0.01, number_of_term)
        report_get_Θ= @report_opt ThermalSpreadingResistance._get_Θ(ThermalSpreadingResistance.Isotropic(), 4.3, 3.2,x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], 0.01, number_of_term)
        sol= ThermalSpreadingResistance.solve(ThermalSpreadingResistance.Isotropic(), x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], number_of_term)
        report_solve = @report_opt sol.Rₛ

        num_issues_Θ_avg = length(JET.get_reports(report_Θ_avg))
        num_issues_get_Θ = length(JET.get_reports(report_get_Θ))
        num_issues_solve = length(JET.get_reports(report_solve))

        @test num_issues_Θ_avg ≈ 0
        @test num_issues_get_Θ ≈ 0
        @test num_issues_solve ≈ 0

        #Compound stability tests
        𝑘ₚ = 5.0
        δₚ = 0.0001
        report_Θ_avg_compound = @report_opt ThermalSpreadingResistance._Θ_avg(Compound(), x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], 0.01, number_of_term)
        report_get_Θ_compound = @report_opt ThermalSpreadingResistance._get_Θ(ThermalSpreadingResistance.Compound(), 4.3, 3.2, x[1], x[2], x[3], x[4], x[5], x[6], x[7], 𝑘ₚ, δₚ, x[8], x[9], x[10],0.01, 100)
        sol_compound = ThermalSpreadingResistance.solve(ThermalSpreadingResistance.Compound(), x[1], x[2], x[3], x[4], x[5], x[6], x[7], 𝑘ₚ, δₚ, x[8], x[9], x[10], number_of_term)
        report_solve_compound = @report_opt sol_compound.Rₛ

        num_issues_Θ_avg_compound = length(JET.get_reports(report_Θ_avg_compound))
        num_issues_get_Θ_compound = length(JET.get_reports(report_get_Θ_compound))
        num_issues_solve_compound = length(JET.get_reports(report_solve_compound))
        
        @test num_issues_Θ_avg_compound ≈ 0
        @test num_issues_get_Θ_compound ≈ 0
        @test num_issues_solve_compound ≈ 0 

end

