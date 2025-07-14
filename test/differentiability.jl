using DifferentiationInterface
using BenchmarkTools
using Mooncake
using ForwardDiff
using Enzyme
using ReverseDiff
include("../src/equations/equations.jl")
f(x::T) where {T<:AbstractVector} = solve(Isotropic(), x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], 50).Rₛ
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
@benchmark ForwardDiff.gradient(f, x)
@benchmark hess = ForwardDiff.hessian(f, x)


Enzyme.gradient(Forward, f, x)
#value_and_gradient(f, AutoForwardDiff(), x) 

value_and_gradient(f, AutoForwardDiff(), x) 
#value_and_gradient(f, AutoEnzyme(),      x)
val, grad = value_and_gradient(f, AutoMooncake(), x)
grad_h = ReverseDiff.gradient(f, x)

# Basic example of using Mooncake to compute gradient and value
using Mooncake

# Define a simple scalar function
g(x) = sum(sin, x) + prod(x)

# Input vector
xg = [1.0, 2.0, 3.0]

# Compute value and gradient using Mooncake
val, grad = value_and_gradient(g, AutoMooncake(), xg)

println("Value of g(x): ", val)
println("Gradient of g(x): ", grad)

# Basic example of using ReverseDiff to compute the gradient

using ReverseDiff

# Define a simple scalar function
h(x) = sum(exp, x) + prod(x)

# Input vector
xh = [1.0, 2.0, 3.0]

# Compute the gradient using ReverseDiff
grad_h = ReverseDiff.gradient(h, xh)

println("Gradient of h(x) using ReverseDiff: ", grad_h)
