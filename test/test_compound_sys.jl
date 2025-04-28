
@testset "CompoundSys" begin

    zeta_val = 1.5
    t1_val = 0.1
    t2_val = 0.2
    k1_val = 10.0
    k2_val = 20.0
    h_val = 100.0
    
    sol = phi_zeta(zeta_val, t1_val, t2_val, k1_val, k2_val, h_val)
    sol_temp = theta_1D_avg(k1_val, k2_val, h_val, t1_val, t2_val)

    @test sol.phi ≈ -2.63605 atol=1e-4
    @test sol_temp.temperature ≈ 3 atol=1e-4
    
end