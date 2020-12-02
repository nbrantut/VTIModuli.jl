using NewtonMethods
using LinearAlgebra: Diagonal, I, Symmetric
using Roots: fzero

include("types.jl")
include("base.jl")
include("invert.jl")
include("cracks.jl")

println("Running test...")

# generate synthetic data, from "true model"
println("Make synthetic data...")

m0 = [0.1,0.05,0.,0.,0.0,deg2rad(28), deg2rad(39), deg2rad(58), deg2rad(90),deg2rad(58),deg2rad(90)]
ρ = 2.5
E₀ = 70
nu₀ = 0.25

C0 = cracks2moduli(E₀, nu₀, m0[1],m0[2],m0[3],m0[4],m0[5])

np=4
nsv=0
nsh=2

args = (np, nsv, nsh)

dobs = _gf(m0, E₀, nu₀, ρ, np, nsv, nsh)

p = [VMeasure(dobs[2k-1],exp(dobs[2k]),.01,0.01) for k in 1:np]
sv = [VMeasure(dobs[2np+2k-1],exp(dobs[2np+2k]),.01,0.01) for k in 1:nsv]
sh = [VMeasure(dobs[2np+2nsv+2k-1],exp(dobs[2np+2nsv+2k]),.01,0.01) for k in 1:nsh]

# a prior model parameters
m1 = [0.0,0.0,0.0,0.0,0.,deg2rad(25), deg2rad(40), deg2rad(55), deg2rad(90),deg2rad(90)]

C1 = cracks2moduli(E₀, nu₀, m1[1],m1[2],m1[3],m1[4],m1[5])

# run inversion
println("Run inversion to find crack parameters...")

m, anglep, anglesv, anglesh, cm = findcracks(p,sv,sh, ρ, E₀, nu₀, (m1[1], 1.0), (m1[2], 1.0), (m1[3], 1e-5), (m1[4], 1e-5), (m1[5],1e-5))

println()
println("Model parameters: input | inverted")
show(IOContext(stdout), "text/plain",[m0[1:5]  m])
println()
println()
dcalc = _gf(vcat(m, anglep, anglesv, anglesh), E₀, nu₀, ρ, np, nsv, nsh)
println("Observations: obs | calc")
show(IOContext(stdout), "text/plain",[dobs  dcalc])
println()

# run inversion moduli
println("Run inversion to find moduli...")

C, anglep, anglesv, anglesh, cm = findmoduli(p,sv,sh, ρ, C1[1,1], C1[1,3], C1[3,3], C1[5,5], C1[6,6], 100)

mc = vcat([C[1,1], C[1,3], C[3,3], C[5,5], C[6,6]], anglep, anglesv, anglesh)

println()
println("Model parameters: input | inverted")
show(IOContext(stdout), "text/plain",C0)
println("\n------")
show(IOContext(stdout), "text/plain",C)
println()
println()
dcalc = _gc(mc, ρ, np, nsv, nsh)
println("Observations: obs | calc")
show(IOContext(stdout), "text/plain",[dobs  dcalc])
