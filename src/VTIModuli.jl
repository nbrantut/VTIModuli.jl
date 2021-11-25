module VTIModuli

using NewtonMethods: quasinewton
using LinearAlgebra: I, Diagonal, Symmetric
using Roots: fzero

export vphaseP, vphaseSV, vphaseSH, moduli2thomsen, velocities0, thomsen2moduli, vgroup, groupangle, dvPdthetav, dvSVdthetav, vgroupSH, groupangleSH, phaseangleP, phaseangleSV, phaseangleSH, vpvsYoungsPoisson, stiffmatrixTI, cracks2moduli
export findthomsen, findmoduli, findcracks
export VMeasure

include("base.jl")
include("types.jl")
include("invert.jl")
include("cracks.jl")

end
