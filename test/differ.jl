
using ForwardDiff
using Enzyme
using ReverseDiff
using JET

function f(x::T) where {T <: AbstractVector}
    solve(Isotropic(), x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], x[9], x[10], 50).Rₛ
end
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
@benchmark Enzyme.gradient(Forward, f, x)
@benchmark ForwardDiff.gradient(f, x)
@benchmark ReverseDiff.gradient(f, x)
@report_opt Enzyme.gradient(Forward, f, x)

@report_opt Enzyme.gradient(Forward, f, x)
@report_opt ForwardDiff.gradient(f, x)
@report_opt ReverseDiff.gradient(f, x)
