using ThermalSpreadingResistance
using Test

@testset "ThermalSpreadingResistance.jl" begin

    include("test_isotropic.jl")
    include("test_compound.jl")
    
end
