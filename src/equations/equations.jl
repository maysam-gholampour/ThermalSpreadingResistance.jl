export Plate, Isotropic, Compound
export solve

abstract type Plate end
struct Isotropic<: Plate end
struct Compound<: Plate end

include("coefficients.jl")
include("theta_avg.jl")
include("isotropic.jl")
include("compound.jl")

