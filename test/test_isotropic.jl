using ThermalSpreadingResistance

a = 0.05
b = 0.05
c = 0.025
d = 0.025
Q = 800.0
𝑘ₛ = 200.0
δₛ = 0.001
hᶜ = 150000.0

using InteractiveUtils: @code_warntype

@code_warntype solve(Isotropic(),a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,15)



@time sol = solve(Isotropic(),a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,100)

sol |> propertynames

sol.Θ_avg
sol.R₁D
sol.Rₛ
sol.Rₜ

sol.Θ(0,0,0)

X = 0:0.001:a
Y = 0:0.001:b
Θ = zeros((length(X),length(Y)))

@inbounds @simd for i in eachindex(X)
    @inbounds @simd for j in eachindex(Y)
        Θ[i,j] = sol.Θ(X[i],Y[j],0.0)
    end
end

using Plots

contourf(X,Y,Θ)

using LinearAlgebra

plot(diag(Θ))