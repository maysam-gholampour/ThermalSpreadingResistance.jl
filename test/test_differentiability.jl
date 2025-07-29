
using ForwardDiff
using Enzyme
using ReverseDiff

@testset "Differentiability" begin
    f(x::T) where {T <:
                   AbstractVector} = solve(
        Isotropic(), x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], 50).Rₛ
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
    enzyme_gradient=Enzyme.gradient(Forward, f, x)
    fowarddiff_gradient=ForwardDiff.gradient(f, x)
    reversediff_gradient = ReverseDiff.gradient(f, x)

    @test enzyme_gradient[1][1] ≈ 16.31377 atol=1e-4
    @test fowarddiff_gradient[1] ≈ 16.31377 atol=1e-4
    @test reversediff_gradient[1] ≈ 16.31377 atol=1e-4
    @test enzyme_gradient[1][3] ≈ -124.7000 atol=1e-4
    @test fowarddiff_gradient[3] ≈ -124.7000 atol=1e-4
    @test reversediff_gradient[3] ≈ -124.7000 atol=1e-4
end
