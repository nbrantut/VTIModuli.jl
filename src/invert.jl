"""
    findthomsen(p::Vector{VMeasure},
                     sv::Vector{VMeasure},
                     sh::Vector{VMeasure},
                     vp₀::Tuple{Real,Real},
                     vs₀::Tuple{Real,Real},
                     ϵ₀::Tuple{Real,Real},
                     δ₀::Tuple{Real,Real},
                     γ₀::Tuple{Real,Real})

    Estimate Thomsen parameters, wavespeeds and phase angles from measurements of group P, SV ans SH wave velocities, using quasinewton method. A priori values ()₀ are tuples with value and standard deviation (a priori gaussian).

Return (vp0, vs0, ϵ, δ, γ), θp, θsv, θsh, CMpost

"""
function findthomsen(p::Vector{VMeasure},
                     sv::Vector{VMeasure},
                     sh::Vector{VMeasure},
                     vp₀::Tuple{Real,Real},
                     vs₀::Tuple{Real,Real},
                     ϵ₀::Tuple{Real,Real},
                     δ₀::Tuple{Real,Real},
                     γ₀::Tuple{Real,Real})
    
    dobs, CDi, np, nsv, nsh = inputdobs(p,sv,sh)
    
    #make vector of model parameters (first guess)
    m₀ = [vp₀[1],vs₀[1],ϵ₀[1],δ₀[1],γ₀[1]]
    σm = [vp₀[2],vs₀[2],ϵ₀[2],δ₀[2],γ₀[2]]
    for x in p
        push!(m₀, phaseangleP(x.angle, vp₀[1]/vs₀[1], ϵ₀[1], δ₀[1]))
        push!(σm, 100.)
    end
    for x in sv
        push!(m₀, phaseangleSV(x.angle, vp₀[1]/vs₀[1], ϵ₀[1], δ₀[1]))
        push!(σm, 100.)
    end
    for x in sh
        push!(m₀, phaseangleSH(x.angle, γ₀[1]))
        push!(σm, 100.)
    end
    
    CM = Diagonal(σm.^2)
    CMi = inv(CM)

    m, CMpost = quasinewton(dobs, _gt, m₀, CMi, CDi, (np,nsv,nsh)) 

    
    return m[1:5], m[6:5+np], m[6+np:5+np+nsv], m[6+np+nsv:5+np+nsv+nsh], CMpost
    
end

"""
    findthomsen(p::Vector{VMeasure},
                     sv::Vector{VMeasure},
                     sh::Vector{VMeasure},
                     vp₀::Real,
                     vs₀::Real,
                     ϵ₀::Real,
                     δ₀::Real,
                     γ₀::Real)

    Method using just numbers as input for a priori values, assuming std=10.0 for all parameters.

"""
function findthomsen(p::Vector{VMeasure},
                     sv::Vector{VMeasure},
                     sh::Vector{VMeasure},
                     vp₀::Real,
                     vs₀::Real,
                     ϵ₀::Real,
                     δ₀::Real,
                     γ₀::Real)

    return findthomsen(p,sv,sh,(vp₀,10.),(vs₀,10.),(ϵ₀,10.),(δ₀,10.),(γ₀,10.))

end

"""
    findmoduli(p::Vector{VMeasure},
                    sv::Vector{VMeasure},
                    sh::Vector{VMeasure},
                    ρ::AbastractFloat,
                    c11₀::Tuple{Real,Real},
                    c13₀::Tuple{Real,Real},
                    c33₀::Tuple{Real,Real},
                    c55₀::Tuple{Real,Real},
                    c66₀::Tuple{Real,Real}

    Estimate stiffness matrix and phase angles from measurements of group P, SV ans SH wave velocities, using quasinewton method. A priori values ()₀ are tuples with value and standard deviation (a priori gaussian).

Return Cijkl, θp, θsv, θsh, CMpost

"""
function findmoduli(p::Vector{VMeasure},
                    sv::Vector{VMeasure},
                    sh::Vector{VMeasure},
                    ρ::Real,
                    c11₀::Tuple{Real,Real},
                    c13₀::Tuple{Real,Real},
                    c33₀::Tuple{Real,Real},
                    c55₀::Tuple{Real,Real},
                    c66₀::Tuple{Real,Real})

    dobs, CDi, np, nsv, nsh = inputdobs(p,sv,sh)
    
    #make vector of model parameters (first guess)
    m₀ = [c11₀[1],c13₀[1],c33₀[1],c55₀[1],c66₀[1]]
    σm = [c11₀[2],c13₀[2],c33₀[2],c55₀[2],c66₀[2]]

    C₀ = stiffmatrixTI(m₀...)
    ϵ₀, δ₀, γ₀ = moduli2thomsen(C₀)
    vp,vs = velocities0(C₀,ρ) # here density does not matter
    
    for x in p
        push!(m₀, phaseangleP(x.angle, vp/vs, ϵ₀, δ₀))
        push!(σm, 100.)
    end
    for x in sv
        push!(m₀, phaseangleSV(x.angle, vp/vs, ϵ₀, δ₀))
        push!(σm, 100.)
    end
    for x in sh
        push!(m₀, phaseangleSH(x.angle, γ₀))
        push!(σm, 100.)
    end
    
    CM = Diagonal(σm.^2)
    CMi = inv(CM)

    m, CMpost = quasinewton(dobs, _gc, m₀, CMi, CDi, (ρ,np,nsv,nsh)) 

    
    return stiffmatrixTI(m[1:5]...), m[6:5+np], m[6+np:5+np+nsv], m[6+np+nsv:5+np+nsv+nsh], CMpost
end

"""
	findmoduli(p::Vector{VMeasure},
                   sv::Vector{VMeasure},
                   sh::Vector{VMeasure},
                   ρ::Real,
                   vp₀::Real,
                   vs₀::Real,
                   errlogC::Real)

Method with isotropic vp,vs as input.
"""
function findmoduli(p::Vector{VMeasure},
                    sv::Vector{VMeasure},
                    sh::Vector{VMeasure},
                    ρ::Real,
                    vp₀::Real,
                    vs₀::Real,
                    errlogC::Real)

    C0 = stiffmatrixTI(vpvsYoungsPoisson(vp₀,vs₀,ρ))

    return findmoduli(p,sv,sh,ρ,
                      (C0[1,1],errlogC),
                      (C0[1,3],errlogC),
                      (C0[3,3],errlogC),
                      (C0[5,5],errlogC),
                      (C0[6,6],errlogC))
    
end
"""
    findmoduli(p::Vector{VMeasure},
                    sv::Vector{VMeasure},
                    sh::Vector{VMeasure},
                    ρ::Real,
                    c11₀::Real,
                    c13₀::Real,
                    c33₀::Real,
                    c55₀::Real,
                    c66₀::Real,
                    errlogC::Real)

Method using just a unique error on moduli.
"""
function findmoduli(p::Vector{VMeasure},
                    sv::Vector{VMeasure},
                    sh::Vector{VMeasure},
                    ρ::Real,
                    c11₀::Real,
                    c13₀::Real,
                    c33₀::Real,
                    c55₀::Real,
                    c66₀::Real,
                    errlogC::Real)

    return findmoduli(p,sv,sh,ρ,
                      (c11₀,errlogC),
                      (c13₀,errlogC),
                      (c33₀,errlogC),
                      (c55₀,errlogC),
                      (c66₀,errlogC))

end

"""
	inputdobs(p,sv,sh)
"""
function inputdobs(p,sv,sh)
    np = length(p)
    nsv = length(sv)
    nsh = length(sh)
    
    # make vector of data and data covariance
    dobs = Float64[]
    σd = Float64[]
    for x in vcat(p,sv,sh)
        append!(dobs, [x.angle, log(x.value)])
        append!(σd, [x.err_angle, x.err_value])
    end
    
    CD = Diagonal(σd.^2)
    CDi = inv(CD)

    return dobs, CDi, np, nsv, nsh
end

function _gc(m,ρ,np,nsv,nsh)
    C₀ = stiffmatrixTI(m[1:5]...)
    ϵ, δ, γ = moduli2thomsen(C₀)
    vp0,vs0 = velocities0(C₀,ρ)
    mt = vcat([vp0,vs0,ϵ,δ,γ], m[6:end])
    return _gt(mt, np, nsv, nsh)
end

function _gt(m,np,nsv,nsh)
    vp0 = m[1]
    vs0 = m[2]
    ϵ = m[3]
    δ = m[4]
    γ = m[5]
    θp = m[6:5+np]
    θsv= m[6+np:5+np+nsv]
    θsh= m[6+np+nsv:5+np+nsv+nsh]

    d = zeros(eltype(m), 2np+2nsv+2nsh)  #NEED eltype HERE to work with ForwardDiff!
    c = 1
    for x in θp
        vp = vphaseP(x,vp0,vs0,ϵ,δ)
        dvvp = dvPdthetav(x,vp0/vs0,ϵ,δ)
        d[c] = groupangle(x,dvvp)
        d[c+1] = log(vgroup(vp, dvvp))
        c+=2
    end
    for x in θsv
        vsv = vphaseSV(x,vp0,vs0,ϵ,δ)
        dvvsv = dvSVdthetav(x,vp0/vs0,ϵ,δ)
        d[c] = groupangle(x,dvvsv)
        d[c+1] = log(vgroup(vsv, dvvsv))
        c+=2
    end
    for x in θsh
        d[c] = groupangleSH(x,γ)
        d[c+1] = log(vgroupSH(d[c], vs0, γ))
        c+=2
    end

    return d
end
