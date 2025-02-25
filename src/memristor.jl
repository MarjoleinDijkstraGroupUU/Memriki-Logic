using DrWatson
@quickactivate

using SpecialFunctions

const ee = 1.6022e-19 # elementary charge
const R_t0 = 1.0 # normalization constant
const R_t = 50e-9
const R_b = 200e-9 / R_t0
const ΔR = R_b - R_t
const qR = R_t / R_b
const RR = R_b * R_t
const L = 10e-6 / R_t0
const T = 293.15 # room temperature
const k_B = 1.3806503e-23 # Boltzman's constant
const kT = k_B * T
const Dc = 1.75e-9 # Diffusion Coefficient of ions in water
const Nₐ = 6.02214076e23 # Avogadros Number
const ρ_b = 0.1e-3 * Nₐ / (1e-3) # bulk salt concentration 0.1mM


const σ = -0.0015 / 1e-18 # wall charge density
const η = 1.01e-3 # dynamic viscosity of water
const ϵ = 0.71e-9 # electric permitivity
const λ_D = sqrt(ϵ * kT / (2ee^2*ρ_b)) # Debyle length
const λ_B = ee^2 / (4π * ϵ * kT) # Bjerrum length

const ψ = 2kT / ee * asinh(2π * λ_D * λ_B * σ) # surface potential
const w = ee * Dc * η / (kT * ϵ * ψ) # ionic mobility / osmotic mobility ratio
const Du = σ / (2ρ_b * R_t) #Dukhin Number

const Δg = -2w * ΔR / R_b * Du
const g₀ = pi * RR / L * 2ρ_b * ee^2 * Dc / (kT) # steady state conductance

const τ = L^2 / 12Dc


Pe(V) = -ee * V * R_b / (kT * R_t * w)
Pe′(V) = -ee * R_b / (kT * R_t * w)

function diffs(h)
"""
test whether numerical derivatives of g∞ coincide with analytical calculation
"""
    println("cdiff : $(cdiff(h))")
    println("cdiff2: $(cdiff2(h))")
    println("bdiff : $(bdiff(h))")
    println("bdiff2: $(bdiff2(h))")
    println("fdiff : $(fdiff(h))")
    println("fdiff2: $(fdiff2(h))")
end

sigmoid(x, tau=1.0) = 1 / (1 + exp(-(x-0.5) / tau)) + 0.0

function g∞(V)
"""
    Calculate g∞ as eq. 5 of https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.130.268401
"""
    P = Pe(V)
    if abs(P) <= 1e-7
        1.0
    else
        Ei = -P * R_t^2 / ΔR * (expinti(P * R_t / ΔR) - expinti(P * qR * R_t / ΔR))
        peclet_part = 1 / (exp(P * qR) - 1) * (exp(-P * qR * R_t / ΔR) / ΔR * (R_t * exp(P * R_t / ΔR) - R_b * exp(P * qR * R_t / ΔR) + Ei) + 1)
        1 + Δg * (R_t / ΔR * (-R_b / ΔR * log(qR) - 1) + peclet_part)
    end
end

function gs(V)
"""
    Calculate g∞ and it's derivative as eq. 5 of https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.130.268401
derivative at 0V is calculated using numerical approximation
"""
    if abs(V) > 7
        println(V)
        return first(gs(7sign(V))),0.0
    end
    P = Pe(V)
    ePR2 = exp(P * R_t / ΔR)
    ePR4 = exp(P * qR * R_t / ΔR)
    Δe = ePR2 - ePR4
    ΔEi = expinti(P * R_t / ΔR) - expinti(P * qR * R_t / ΔR)
    peclet_part = 1 / (exp(P * qR) - 1) * (1/(ePR4*ΔR)*(R_t * ePR2 - R_b * ePR4 - P * R_t^2 / ΔR * ΔEi) + 1)
    if P == 0.0
        g = 1.0
        g′ = -1.02426
    else
        g = 1.0 + Δg * (R_t / ΔR * (-R_b / ΔR * log(qR) - 1) + peclet_part)
        # g′ = Pe′(V) * Δg * R_t / (ΔR^2 * Δe) * ((qR * ePR4 - ePR2) * (peclet_part) - R_t * ΔEi) # has divergence problems when V -> 0
        g′ = (g∞(V+1e-5)-g∞(V-1e-5))/2e-5
    end
    (g, g′)
end
