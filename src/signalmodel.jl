"""
    stfr(M0, T1, T2, Δω, κ, Tfree, Tg, α, β, ϕ, TE, [spoil])

Compute the steady-state signal generated from an STFR scan of the given tissue
parameters with the given scan parameters, assuming instantaneous RF pulses and
ideal spoiling.

# Arguments
- `M0::Real`: Equilibrium magnetization
- `T1::Real`: Spin-lattice recovery time constant (ms)
- `T2::Real`: Spin-spin recovery time constant (ms)
- `Δω::Real`: Off-resonance frequency (rad/s)
- `κ::Real`: Variation in tip angles
- `Tfree::Real`: Free precession time (ms)
- `Tg::Real`: Spoiler gradient time (ms)
- `α::Real`: Tip-down angle (rad)
- `β::Real`: Tip-up angle (rad)
- `ϕ::Real`: Phase of tip-up RF pulse (rad)
- `TE::Real`: Echo time (ms)
- `spoil::Val = Val(true)`: Whether to include RF spoiling

# Return
- `signal::Complex`: Signal generated from the STFR scan
"""
function stfr(
    M0::Real,
    T1::Real,
    T2::Real,
    Δω::Real,
    κ::Real,
    Tfree::Real,
    Tg::Real,
    α::Real,
    β::Real,
    ϕ::Real,
    TE::Real,
    ::Val{spoil} = Val(true)
) where {spoil}

    if spoil
        return   stfr_spoil(M0, T1, T2, Δω, κ, Tfree, Tg, α, β, ϕ, TE)
    else
        return stfr_unspoil(M0, T1, T2, Δω, κ, Tfree, Tg, α, β, ϕ, TE)
    end

end

"""
    stfr(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ, Tfree, Tg, α, β, ϕ, TE, [spoil])

Compute the steady-state signal generated from an STFR scan of the given
two-compartment tissue parameters with the given scan parameters, assuming
instantaneous RF pulses, ideal spoiling, and no exchange between the two
compartments.

# Arguments
- `M0::Real`: Equilibrium magnetization
- `ff::Real`: Fraction of tissue that is fast compartment
- `T1f::Real`: Fast compartment spin-lattice recovery time constant (ms)
- `T1s::Real`: Slow compartment spin-lattice recovery time constant (ms)
- `T2f::Real`: Fast compartment spin-spin recovery time constant (ms)
- `T2s::Real`: Slow compartment spin-spin recovery time constant (ms)
- `Δωf::Real`: Fast compartment off-resonance frequency (added to Δω) (rad/s)
- `Δω::Real`: Bulk off-resonance frequency (rad/s)
- `κ::Real`: Variation in tip angles
- `Tfree::Real`: Free precession time (ms)
- `Tg::Real`: Spoiler gradient time (ms)
- `α::Real`: Tip-down angle (rad)
- `β::Real`: Tip-up angle (rad)
- `ϕ::Real`: Phase of tip-up RF pulse (rad)
- `TE::Real`: Echo time (ms)
- `spoil::Val = Val(true)`: Whether to include RF spoiling

# Return
- `signal::Complex`: Signal generated from the two-compartment STFR scan
"""
function stfr(
    M0::Real,
    ff::Real,
    T1f::Real,
    T1s::Real,
    T2f::Real,
    T2s::Real,
    Δωf::Real,
    Δω::Real,
    κ::Real,
    Tfree::Real,
    Tg::Real,
    α::Real,
    β::Real,
    ϕ::Real,
    TE::Real,
    ::Val{spoil} = Val(true)
) where {spoil}

    if spoil
        return   ff  * stfr_spoil(M0, T1f, T2f, Δω + Δωf,
                                  κ, Tfree, Tg, α, β, ϕ, TE) +
            (1 - ff) * stfr_spoil(M0, T1s, T2s, Δω,
                                  κ, Tfree, Tg, α, β, ϕ, TE)
    else
        return ff  * stfr_unspoil(M0, T1f, T2f, Δω + Δωf,
                                  κ, Tfree, Tg, α, β, ϕ, TE) +
          (1 - ff) * stfr_unspoil(M0, T1s, T2s, Δω,
                                  κ, Tfree, Tg, α, β, ϕ, TE)
    end

end

"Compute the STFR signal using the spoiled STFR signal model."
function stfr_spoil(M0, T1, T2, Δω, κ, Tfree, Tg, α, β, ϕ, TE)

    κα = κ * α
    κβ = κ * β
    θfree = Δω * Tfree/1000

    n = exp(-Tg/T1) * (1 - exp(-Tfree/T1)) * sin(κα) * cos(κβ) +
        (1 - exp(-Tg/T1)) * sin(κα)

    d = 1 - exp(-Tg/T1 - Tfree/T2) * sin(κα) * sin(κβ) * cos(θfree - ϕ) -
            exp(-Tg/T1 - Tfree/T1) * cos(κα) * cos(κβ)

    return (M0 * n / d) * exp(-TE/T2 - im*Δω*TE/1000)

end

"Compute the STFR signal using the unspoiled STFR signal model."
function stfr_unspoil(M0, T1, T2, Δω, κ, Tfree, Tg, α, β, ϕ, TE)

    κα = κ * α
    κβ = κ * β
    θfree = Δω * Tfree/1000
    Ef1 = exp(-Tfree/T1)
    Ef2 = exp(-Tfree/T2)
    Eg1 = exp(-Tg/T1)
    Eg2 = exp(-Tg/T2)

    a = -im * Eg2 * (Ef2 * (-1 + Eg1 + (-1 + Ef1) * Eg1 * cos(κβ)) *
        cos(θfree - ϕ) * sin(κα) + (Ef1 * (-1 + Eg1) + (-1 + Ef1) * cos(κα)) *
        sin(κβ) + im * Ef2 * (-1 + Eg1 + (-1 + Ef1) * Eg1 * cos(κβ)) * sin(κα) *
        sin(θfree - ϕ))

    b = Eg2 * (Ef2 * ((-1 + Ef1) * Eg1 + (-1 + Eg1) * cos(κβ)) * cos(θfree - ϕ) *
        sin(κα) - (-1 + Ef1 + Ef1 * (-1 + Eg1) * cos(κα)) * sin(κβ) + im * Ef2 *
        ((-1 + Ef1) * Eg1 + (-1 + Eg1) * cos(κβ)) * sin(κα) * sin(θfree - ϕ))

    c = im * ((-1 + Eg1 + (-1 + Ef1) * Eg1 * cos(κβ)) * sin(κα) + Ef2 * Eg2^2 *
        (Ef1 * (-1 + Eg1) + (-1 + Ef1) * cos(κα)) * sin(κβ) * (cos(θfree - ϕ) +
        im * sin(θfree - ϕ)))

    d = Eg2 * (-Ef2 * (-1 + Ef1 * Eg1) * (1 + cos(κα) * cos(κβ)) * cos(θfree - ϕ)
        + (Ef1 - Ef2^2 * Eg1) * sin(κα) * sin(κβ))

    e = Ef2 * (-1 + Ef1 * Eg1) * Eg2 * (cos(κα) + cos(κβ)) * sin(θfree - ϕ)

    f = -1 + Ef1 * Ef2^2 * Eg1 * Eg2^2 + (Ef1 * Eg1 - Ef2^2 * Eg2^2) * cos(κα) *
        cos(κβ) + Ef2 * (Eg1 - Ef1 * Eg2^2) * cos(θfree - ϕ) * sin(κα) * sin(κβ)

    s = sqrt(f^2 - d^2 - e^2)

    tmp1 = -c / s
    tmp2 = (a * d + b * e) / (d^2 + e^2)
    tmp3 = (-f - s) / s

    return M0 * (tmp1 - tmp2 * tmp3) * exp(-TE/T2 - im*Δω*TE/1000)

end
