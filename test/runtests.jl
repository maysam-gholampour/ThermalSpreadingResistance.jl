using ThermalSpreadingResistance
using Test

@testset "ThermalSpreadingResistance.jl" begin
    include("test_isotropic.jl")
    include("test_compound.jl")
    include("test_differentiability.jl")
    include("test_stability.jl")
end
