using NewtonMethods
using LinearAlgebra: Diagonal, I, Symmetric
using Roots: fzero

include("types.jl")
include("base.jl")
include("invert.jl")

println("Running test...")
println("Make synthetic data...")

m0 = [2,1.4,0.1,0.15,0.03,deg2rad(25), deg2rad(40), deg2rad(55), deg2rad(90),deg2rad(90)]
ρ = 1.

C0 = thomsen2moduli(m0[1],m0[2],ρ,m0[3],m0[4],m0[5])

np=4
nsv=0
nsh=1

args = (np, nsv, nsh)

dobs = _gt(m0, np, nsv, nsh)

p = [VMeasure(dobs[2k-1],exp(dobs[2k]),.1,0.1) for k in 1:np]
sh = [VMeasure(dobs[2np+2k-1],exp(dobs[2np+2k]),.1,0.1) for k in 1:nsh]

# a prior moel parameters
m1 = [2.3,1.7,0.0,0.0,0.,deg2rad(25), deg2rad(40), deg2rad(55), deg2rad(90),deg2rad(90)]

C1 = thomsen2moduli(m1[1],m1[2],ρ,m1[3],m1[4],m1[5])

# run inversion thomsen
println("Run inversion to find Thomsen parameters...")

m, anglep, anglesv, anglesh, cm = findthomsen(p,VMeasure[],sh, m1[1], m1[2], m1[3], m1[4], m1[5])

println()
println("Model parameters: input | inverted")
show(IOContext(stdout), "text/plain",[m0[1:5]  m])
println()
println()
dcalc = _gt(vcat(m, anglep, anglesv, anglesh), np, nsv, nsh)
println("Observations: obs | calc")
show(IOContext(stdout), "text/plain",[dobs  dcalc])
println()

# run inversion moduli
println("Run inversion to find moduli...")

C, anglep, anglesv, anglesh, cm = findmoduli(p,VMeasure[],sh, ρ, C1[1,1], C1[1,3], C1[3,3], C1[5,5], C1[6,6], 100)

println()
println("Model parameters: input | inverted")
show(IOContext(stdout), "text/plain",C0)
println("\n------")
show(IOContext(stdout), "text/plain",C)
println()
println()
dcalc = _gc(vcat([C[1,1], C[1,3], C[3,3], C[5,5], C[6,6]], anglep, anglesv, anglesh), ρ, np, nsv, nsh)
println("Observations: obs | calc")
show(IOContext(stdout), "text/plain",[dobs  dcalc])


#=
function testgradient(x)
    r = zeros(eltype(x),2length(x))
    for (c,elem) in enumerate(x)
        r[2c-1] = elem
        r[2c] = atan(exp(elem))
    end
    return r
end
=#
