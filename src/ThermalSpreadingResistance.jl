module ThermalSpreadingResistance

    using PrecompileTools: @setup_workload, @compile_workload  

    include("equations/equations.jl")

    @setup_workload begin
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

        @compile_workload begin
            sol = solve(Isotropic(),a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ,100)

            sol = solve(Compound(),a,b,c,d,Q,𝑘ₛ,δₛ,𝑘ₚ,δₚ,hᶜ,Xᶜ,Yᶜ,100)
        end
    end

end
