"""
	cracks2moduli(E₀, nu₀, α1, α3, β11, β13, β33)

Compute 6×6 stiffness matrix (in Voight notation) for a cracked material characterised by normalised crack density parameters α1, α3, β11, β13 and β33, with solid matrix moduli E₀ and nu₀. Crack density parameters are normalised by h = 3E₀(2-nu₀)/(32(1-nu₀²)). From Sayers and Kachanov 1995.
"""
function cracks2moduli(E₀, nu₀, α1, α3, β11, β13, β33)
    S011 = 1/E₀
    S012 = -nu₀/E₀
    h = cracknormfactor(E₀, nu₀)

    D = (S011 + α3/h + β33/h)*(S011 + S012 + α1/h + 4β11/h/3) -
        2*(S012 + β13/h)^2

    c11 = (1/2)*((S011+α3/h+β33/h)/D + 1/(S011-S012+α1/h+2β11/h/3))
    c33 = (S011 + S012 + α1/h + 4β11/h/3)/D
    c44 = 1/(2S011 - 2S012 + α1/h + α3/h + 4β13/h)
    c13 = -(S012 + β13/h)/D
    c66 = 1/(S011-S012+α1/h+2β11/h/3)/2

    return stiffmatrixTI(c11, c13, c33, c44, c66)
end

"""
	cracknormfactor(E₀, ν₀)

Compute normalisation factor h=3E₀*(2-nu₀)/(32(1-nu₀²)).
"""
function cracknormfactor(E₀, ν₀)
    return 3E₀*(2-ν₀)/(32(1-ν₀^2))
end

"""
	findcracks(p::Vector{VMeasure},
                    sv::Vector{VMeasure},
                    sh::Vector{VMeasure},
                    ρ::Real,
                    E₀::Real,
                    nu₀::Real,
                    α₁₀::Tuple{Real,Real},
                    α₃₀::Tuple{Real,Real},
                    β₁₁₀::Tuple{Real,Real},
                    β₁₃₀::Tuple{Real,Real},
                    β₃₃₀::Tuple{Real,Real})

Estimate normalised crack density parameters and phase angles from measurements of group P, SV ans SH wave velocities, using quasinewton method. A priori values ()₀ are tuples with value and standard deviation (a priori gaussian).

Return α₁₁, α₃₃, β₁₁₁₁, β₁₁₃₃, β₃₃₃₃, θp, θsv, θsh, CMpost.
"""
function findcracks(p::Vector{VMeasure},
                    sv::Vector{VMeasure},
                    sh::Vector{VMeasure},
                    ρ::Real,
                    E₀::Real,
                    nu₀::Real,
                    α₁₀::Tuple{Real,Real},
                    α₃₀::Tuple{Real,Real},
                    β₁₁₀::Tuple{Real,Real},
                    β₁₃₀::Tuple{Real,Real},
                    β₃₃₀::Tuple{Real,Real})

    dobs, CDi, np, nsv, nsh = inputdobs(p,sv,sh)
    
    #make vector of model parameters (first guess)
    m₀ = [α₁₀[1],α₃₀[1],β₁₁₀[1],β₁₃₀[1],β₃₃₀[1]]
    σm = [α₁₀[2],α₃₀[2],β₁₁₀[2],β₁₃₀[2],β₃₃₀[2]]

    C₀ = cracks2moduli(E₀, nu₀, m₀...)
    ϵ₀, δ₀, γ₀ = moduli2thomsen(C₀)
    vp,vs = velocities0(C₀,ρ) # here density does not matter
    
    for x in p
        push!(m₀, phaseangleP(x.angle, vp/vs, ϵ₀, δ₀))
        push!(σm, 100)
    end
    for x in sv
        push!(m₀, phaseangleSV(x.angle, vp/vs, ϵ₀, δ₀))
        push!(σm, 100)
    end
    for x in sh
        push!(m₀, phaseangleSH(x.angle, γ₀))
        push!(σm, 100)
    end
    
    CM = Diagonal(σm.^2)
    CMi = inv(CM)

    m, CMpost = quasinewton(dobs, _gf, m₀, CMi, CDi, (E₀, nu₀, ρ,np,nsv,nsh)) 

    
    return m[1:5], m[6:5+np], m[6+np:5+np+nsv], m[6+np+nsv:5+np+nsv+nsh], CMpost
    
end

function _gf(m,E₀,nu₀,ρ,np,nsv,nsh)
    C₀ = cracks2moduli(E₀, nu₀, m[1:5]...)
    ϵ, δ, γ = moduli2thomsen(C₀)
    vp0,vs0 = velocities0(C₀,ρ)
    mt = vcat([vp0,vs0,ϵ,δ,γ], m[6:end])
    return _gt(mt, np, nsv, nsh)
end
