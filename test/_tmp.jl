
using ThermalSpreadingResistance
using ThermalSpreadingResistance: _Θ_avgr

a = 0.05
b = 0.05
c = 0.025
d = 0.025
Q = 800.0
𝑘ₛ = 200.0
δₛ = 0.001
hᶜ = 150000.0

Xᶜ = 0.5 * a
Yᶜ = 0.5 * b

number_of_term = 50



Θ₁D = 0.01

#function _Θ_avgr(a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ,Θ₁D,number_of_term)

_Θ_avgr(a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ,0.01,number_of_term)

@code_warntype _Θ_avgr(a,b,c,d,Q,𝑘ₛ,δₛ,hᶜ,Xᶜ,Yᶜ,0.01,number_of_term)

using BenchmarkTools

@benchmark _Θ_avgr($a,$b,$c,$d,$Q,$𝑘ₛ,$δₛ,$hᶜ,$Xᶜ,$Yᶜ,0.01,$number_of_term)