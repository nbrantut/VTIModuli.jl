function vphaseP(θ,VP0,VS0,ϵ,δ)
    f = 1 - (VS0^2)/(VP0^2)
    D = sqrt(1 + 4*sin(θ)^2/f*(2*δ*cos(θ)^2 - ϵ*cos(2θ)) +
             (2*ϵ*sin(θ)^2/f)^2)
    vsq = 1 + ϵ*sin(θ)^2 - f/2 + D*f/2
    return sqrt(vsq)*VP0
end

function vphaseSV(θ,VP0,VS0,ϵ,δ)
    f = 1 - (VS0^2)/(VP0^2)
    D = sqrt(1 + 4*sin(θ)^2/f*(2*δ*cos(θ)^2 - ϵ*cos(2θ)) +
             (2*ϵ*sin(θ)^2/f)^2)
    vsq = 1 + ϵ*sin(θ)^2 - f/2 - D*f/2
    return sqrt(vsq)*VP0
end

function vphaseSH(θ,VS0,γ)
    return VS0*sqrt(1 + 2*γ*sin(θ)^2)
end

function vphaseP(θ,c,ρ)
    stsq = sin(θ)^2
    ctsq = cos(θ)^2
    D = sqrt(((c[1,1]-c[5,5])*stsq - (c[3,3]-c[5,5])*ctsq)^2 +
             4*(c[1,3]+c[5,5])^2*stsq*ctsq)
    rv = (c[1,1] + c[5,5])*stsq + (c[3,3]+c[5,5])*ctsq + D
    return sqrt(rv/(2ρ))
end

function vphaseSV(θ,c,ρ)
    stsq = sin(θ)^2
    ctsq = cos(θ)^2
    D = sqrt(((c[1,1]-c[5,5])*stsq - (c[3,3]-c[5,5])*ctsq)^2 +
             4*(c[1,3]+c[5,5])^2*stsq*ctsq)
    rv = (c[1,1] + c[5,5])*stsq + (c[3,3]+c[5,5])*ctsq - D
    return sqrt(rv/(2ρ))
end

function vphaseSH(θ,c,ρ)
    return sqrt((c[6,6]*sin(θ)^2 + c[5,5]*cos(θ)^2)/ρ)
end

function moduli2thomsen(c)
    epsilon = (c[1,1] - c[3,3])/(2*c[3,3])
    delta = ((c[1,3]+c[5,5])^2 - (c[3,3]-c[5,5])^2)/(2*c[3,3]*(c[3,3]-c[5,5]))
    gamma = (c[6,6]-c[5,5])/(2*c[5,5])
    return epsilon, delta, gamma
end

function velocities0(c,ρ)
    return (sqrt(c[3,3]/ρ), sqrt(c[5,5]/ρ))
end

function thomsen2moduli(VP0,VS0,ρ,ϵ,δ,γ)
    c33 = ρ*VP0^2
    c55 = ρ*VS0^2
    c11 = c33*(2*ϵ+1)
    c66 = c55*(2*γ+1)
    c13 = sqrt(2*c33*δ*(c33-c55) + (c33-c55)^2) - c55

    return stiffmatrixTI(c11,c13,c33,c55,c66)
end

function vpvsYoungsPoisson(vp,vs,ρ)

    nu = (0.5*(vp/vs)^2 - 1)/((vp/vs)^2 - 1)
    E = ρ*vs^2*2*(1+nu)

    return E, nu
end

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
    return Symmetric(c)
end


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

function dvPdthetav(θ,VP0VS0,ϵ,δ)

    f = 1 - VP0VS0^(-2)

    a = sqrt((δ+f-δ*cos(4θ)-4*ϵ*cos(2θ)*sin(θ)^2+4*ϵ^2*sin(θ)^4)/f)
    
    vpv = (δ*sin(4θ) + ϵ*sin(2θ)*(1+ϵ-(2+ϵ)*cos(2θ) + a))/(a*(2+2*ϵ*sin(θ)^2 + f*(-1+a)))

    return vpv
end

function dvSVdthetav(θ,VP0VS0,ϵ,δ)

    f = 1 - VP0VS0^(-2)

    a = sqrt((δ+f-δ*cos(4θ)-4*ϵ*cos(2θ)*sin(θ)^2+4*ϵ^2*sin(θ)^4)/f)
    
    vpv = (δ*sin(4θ) - ϵ*sin(2θ)*(-1-ϵ+(2+ϵ)*cos(2θ) + a))/(a*(-2-2*ϵ*sin(θ)^2 + f*(1+a)))

    return vpv
end

function vgroup(v,dvv)
    return v*sqrt(1 + dvv^2)
end

function groupangle(θ,dvdθv)
    return atan((sin(θ) + dvdθv*cos(θ))/(cos(θ) - dvdθv*sin(θ)))
end

function vgroupSH(ψ,VS0,γ)
    return VS0*sqrt((1+2*γ)/(1+2*γ*cos(ψ)^2))
end

function groupangleSH(θ,γ)
    return atan(tan(θ)*(1+2*γ))
end

function phaseangleP(ψ,VP0VS0,ϵ,δ)
    fzero(θ->(ψ - groupangle(θ,dvPdthetav(θ,VP0VS0,ϵ,δ))),ψ)
end

function phaseangleSV(ψ,VP0VS0,ϵ,δ)
    fzero(θ->(ψ - groupangle(θ,dvSVdthetav(θ,VP0VS0,ϵ,δ))),ψ)
end

function phaseangleSH(ψ,γ)
    return atan(tan(ψ)/(1+2*γ))
end
