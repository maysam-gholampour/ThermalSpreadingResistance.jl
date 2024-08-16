using ThermalSpreadingResistance

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

using InteractiveUtils: @code_warntype

@code_warntype solve(Compound(),a,b,c,d,Q,𝑘ₛ,δₛ,𝑘ₚ,δₚ,hᶜ,15)



@time sol = solve(Compound(),a,b,c,d,Q,𝑘ₛ,δₛ,𝑘ₚ,δₚ,hᶜ,100)

sol |> propertynames

sol.Θ_avg
sol.R₁D
sol.Rₛ
sol.Rₜ
