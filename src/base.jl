"""
    vphaseP(θ, VP0, VS0, ϵ, δ)

Compute P wave phase velocity along phase angle θ (with respect to axis of symmetry) based on Vp₀ and Vs₀ (P and S wave speed along axis of symmetry) and Thomsen parameters ϵ and δ.
"""
function vphaseP(θ,VP0,VS0,ϵ,δ)
    f = 1 - (VS0^2)/(VP0^2)
    D = sqrt(1 + 4*sin(θ)^2/f*(2*δ*cos(θ)^2 - ϵ*cos(2θ)) +
             (2*ϵ*sin(θ)^2/f)^2)
    vsq = 1 + ϵ*sin(θ)^2 - f/2 + D*f/2
    return sqrt(vsq)*VP0
end

"""
    vphaseSV(θ, VP0, VS0, ϵ, δ)

Compute SV wave phase velocity along phase angle θ (with respect to axis of symmetry) based on Vp₀ and Vs₀ (P and S wave speed along axis of symmetry) and Thomsen parameters ϵ and δ.
"""
function vphaseSV(θ,VP0,VS0,ϵ,δ)
    f = 1 - (VS0^2)/(VP0^2)
    D = sqrt(1 + 4*sin(θ)^2/f*(2*δ*cos(θ)^2 - ϵ*cos(2θ)) +
             (2*ϵ*sin(θ)^2/f)^2)
    vsq = 1 + ϵ*sin(θ)^2 - f/2 - D*f/2
    return sqrt(vsq)*VP0
end

"""
    vphaseSH(θ, VP0, VS0, γ)

Compute SH wave phase velocity along phase angle θ (with respect to axis of symmetry) based on Vp₀ and Vs₀ (P and S wave speed along axis of symmetry) and Thomsen parameter γ.
"""
function vphaseSH(θ,VS0,γ)
    return VS0*sqrt(1 + 2*γ*sin(θ)^2)
end

"""
    vphaseP(θ, c::AbstractMatrix, ρ)

Compute P wave phase velocity along phase angle θ (with respect to axis of symmetry) based on stiffness matrix c (in Voight notation) and density ρ.
"""
function vphaseP(θ,c::AbstractMatrix,ρ)
    stsq = sin(θ)^2
    ctsq = cos(θ)^2
    D = sqrt(((c[1,1]-c[5,5])*stsq - (c[3,3]-c[5,5])*ctsq)^2 +
             4*(c[1,3]+c[5,5])^2*stsq*ctsq)
    rv = (c[1,1] + c[5,5])*stsq + (c[3,3]+c[5,5])*ctsq + D
    return sqrt(rv/(2ρ))
end

"""
    vphaseSV(θ, c::AbstractMatrix, ρ)

Compute SV wave phase velocity along phase angle θ (with respect to axis of symmetry) based on stiffness matrix c (in Voight notation) and density ρ.
"""
function vphaseSV(θ,c::AbstractMatrix,ρ)
    stsq = sin(θ)^2
    ctsq = cos(θ)^2
    D = sqrt(((c[1,1]-c[5,5])*stsq - (c[3,3]-c[5,5])*ctsq)^2 +
             4*(c[1,3]+c[5,5])^2*stsq*ctsq)
    rv = (c[1,1] + c[5,5])*stsq + (c[3,3]+c[5,5])*ctsq - D
    return sqrt(rv/(2ρ))
end

"""
    vphaseSH(θ, c::AbstractMatrix, ρ)

Compute SH wave phase velocity along phase angle θ (with respect to axis of symmetry) based on stiffness matrix c (in Voight notation) and density ρ.
"""
function vphaseSH(θ,c::AbstractMatrix,ρ)
    return sqrt((c[6,6]*sin(θ)^2 + c[5,5]*cos(θ)^2)/ρ)
end

"""
    moduli2thomsen(c::AbstractMatrix)

Compute Thomsen's parameters (ϵ, δ, γ) based on stiffness matrix c (in Voight notation).
"""
function moduli2thomsen(c::AbstractMatrix)
    epsilon = (c[1,1] - c[3,3])/(2*c[3,3])
    delta = ((c[1,3]+c[5,5])^2 - (c[3,3]-c[5,5])^2)/(2*c[3,3]*(c[3,3]-c[5,5]))
    gamma = (c[6,6]-c[5,5])/(2*c[5,5])
    return epsilon, delta, gamma
end

"""
    velocities0(c::AbstractMatrix, ρ)

Compute phase velocities Vp₀ and Vs₀ (P and S wave speed along symmetry axis) based on stiffness matrix c (in Voight notation) and density ρ.
"""
function velocities0(c::AbstractMatrix,ρ)
    return (sqrt(c[3,3]/ρ), sqrt(c[5,5]/ρ))
end

"""
    thomsen2moduli(VP0, VS0, ρ, ϵ, δ, γ)

Compute 6×6 stiffness matrix from Vp₀, Vs₀ (P and S wave speed along symmetry axis), density ρ, and Thomsen parameters ϵ, δ and γ.
"""
function thomsen2moduli(VP0,VS0,ρ,ϵ,δ,γ)
    c33 = ρ*VP0^2
    c55 = ρ*VS0^2
    c11 = c33*(2*ϵ+1)
    c66 = c55*(2*γ+1)
    c13 = sqrt(2*c33*δ*(c33-c55) + (c33-c55)^2) - c55

    return stiffmatrixTI(c11,c13,c33,c55,c66)
end

"""
    vpvsYoungsPoisson(vp,vs,ρ)

Compute (E, ν), Youngs modulus and Poisson's ratio, from P and S wave speeds, assuming isotropy.
"""
function vpvsYoungsPoisson(vp,vs,ρ)

    nu = (0.5*(vp/vs)^2 - 1)/((vp/vs)^2 - 1)
    E = ρ*vs^2*2*(1+nu)

    return E, nu
end

"""
    stiffmatrixTI(E, nu)

Compute 6×6 stiffness matrix in an isotropic medium characterised by Young's modulus E and Poisson's ratio nu.
"""
function stiffmatrixTI(E,nu)
    c0 = E/(1+nu)/(1-2nu)
    c = zeros(typeof(E),6,6)
    c[1,1] = 1-nu
    c[2,2] = 1-nu
    c[3,3] = 1-nu
    c[1,2] = nu
    c[1,3] = nu
    c[2,3] = nu
    c[4,4] = (1-2nu)/2
    c[5,5] = (1-2nu)/2
    c[6,6] = (1-2nu)/2
    return Symmetric(c).*c0
end

"""
    stiffmatrixTI(c11, c13, c33, c55, c66)

Compute 6×6 stiffness matrix in a VTI medium characterised by components (in Voight notation) c11, c13, c33, c55 and c66.
"""
function stiffmatrixTI(c11,c13,c33,c55,c66)
    c = zeros(typeof(c11),6,6)
    c[1,1] = c11
    c[1,2] = c11-2*c66
    c[1,3] = c13
    c[2,2] = c11
    c[2,3] = c13
    c[3,3] = c33
    c[4,4] = c55
    c[5,5] = c55
    c[6,6] = c66
    return Symmetric(c)
end

"""
    dvPdthetav(θ, VP0VS0, ϵ, δ)

Compute (1/v)∂v/∂θ, derivative of P wave phase velocity  with respect to phase angle θ, normalised by Vp, based on ratio Vp₀/Vs₀ (P and S wave speed along symmetry axis) and Thomsen parameters ϵ and δ.
"""
function dvPdthetav(θ,VP0VS0,ϵ,δ)

    f = 1 - VP0VS0^(-2)

    a = sqrt((δ+f-δ*cos(4θ)-4*ϵ*cos(2θ)*sin(θ)^2+4*ϵ^2*sin(θ)^4)/f)
    
    vpv = (δ*sin(4θ) + ϵ*sin(2θ)*(1+ϵ-(2+ϵ)*cos(2θ) + a))/(a*(2+2*ϵ*sin(θ)^2 + f*(-1+a)))

    return vpv
end

"""
    dvPdthetav(θ, VP0VS0, ϵ, δ)

Compute (1/v)∂v/∂θ, derivative of SV wave phase velocity  with respect to phase angle θ, normalised by Vsv, based on ratio Vp₀/Vs₀ (P and S wave speed along symmetry axis) and Thomsen parameters ϵ and δ.
"""
function dvSVdthetav(θ,VP0VS0,ϵ,δ)

    f = 1 - VP0VS0^(-2)

    a = sqrt((δ+f-δ*cos(4θ)-4*ϵ*cos(2θ)*sin(θ)^2+4*ϵ^2*sin(θ)^4)/f)
    
    vpv = (δ*sin(4θ) - ϵ*sin(2θ)*(-1-ϵ+(2+ϵ)*cos(2θ) + a))/(a*(-2-2*ϵ*sin(θ)^2 + f*(1+a)))

    return vpv
end

"""
    vgroup(v,dvv)

Compute group wave velocity based on phase velocity v and normalised derivative of phase velocity with respect to phase angle dvv.
"""
function vgroup(v,dvv)
    return v*sqrt(1 + dvv^2)
end

"""
    groupangle(θ, dvdθv)

Compute group angle from phase angle θ and normalised derivative of phase velocity with respect to phase angle.
"""
function groupangle(θ,dvdθv)
    return atan((sin(θ) + dvdθv*cos(θ))/(cos(θ) - dvdθv*sin(θ)))
end

"""
    vgroupSH(ψ, VS0, γ)

Compute SH group velocity along group angle ψ, based on Vs₀ (S wave speed along symmetry axis) and Thomsen parameter γ.
"""
function vgroupSH(ψ,VS0,γ)
    return VS0*sqrt((1+2*γ)/(1+2*γ*cos(ψ)^2))
end

"""
    groupangleSH(θ,γ)

Compute SH group angle based on phase angle θ and Thomsen parameter γ.
"""
function groupangleSH(θ,γ)
    return atan(tan(θ)*(1+2*γ))
end

"""
    phaseangleP(ψ, VP0VS0, ϵ, δ)

Determine (using root finding method) the P wave phase angle based on group angle ψ, ratio Vp₀/Vs₀ (P and S wave speed along symmetry axis), and Thomsen parameters ϵ and δ.
"""
function phaseangleP(ψ,VP0VS0,ϵ,δ)
    fzero(θ->(ψ - groupangle(θ,dvPdthetav(θ,VP0VS0,ϵ,δ))),ψ)
end

"""
    phaseangleSV(ψ, VP0VS0, ϵ, δ)

Determine (using root finding method) the SV wave phase angle based on group angle ψ, ratio Vp₀/Vs₀ (P and S wave speed along symmetry axis), and Thomsen parameters ϵ and δ.
"""
function phaseangleSV(ψ,VP0VS0,ϵ,δ)
    fzero(θ->(ψ - groupangle(θ,dvSVdthetav(θ,VP0VS0,ϵ,δ))),ψ)
end

"""
    phaseangleSH(ψ, γ)

Compute SH wave phase angle based on group angle ψ and Thomsen parameter γ.
"""
function phaseangleSH(ψ,γ)
    return atan(tan(ψ)/(1+2*γ))
end
