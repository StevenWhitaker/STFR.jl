"""
    stfrblochsim(M0, T1, T2, Δω, κ, Tfree, Tg, α, β, ϕ[, TE]; how, kwargs...)

Calculate the STFR signal of a single-compartment spin using Bloch simulation.

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
- `TE::Real = Tfree / 2`: Echo time (ms)
- `how::Symbol = :ideal`: Specify what assumptions to make in the Bloch
    simulation; options:
    - `:ideal`: Assume instantaneous RF pulses and ideal spoiling; this option
        accepts no additional keyword arguments
    - `:grad`: Assume instantaneous RF pulses and nonideal gradient spoiling
        (no RF spoiling); this option can accept the following additional
        keyword arguments:
        - `nspins::Integer = 10`: Number of spins to simulate; the returned
            value will be the complex mean of the spin ensemble
        - `ncycles::Real = 1`: Number of cycles of phase caused by the spoiler
            gradient across the `nspins` spins (typically integer-valued)
    - `:rfspoil`: Assume instantaneous RF pulses and nonideal gradient and RF
        spoiling; this option can accept the following additional keyword
        arguments:
        - `nspins` and `ncycles`, which are described under option `:grad`
        - `Δθinc::Real = deg2rad(117)`: Quadratic RF phase increment (rad)
        - `nTR::Integer = 100`: Number of TR's to simulate to produce a pseudo
            steady-state

# Return
- `M::Complex{Float64}`: Signal generated from an STFR scan
"""
function stfrblochsim(
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
    TE::Real = Tfree / 2;
    how::Symbol = :ideal,
    kwargs...
)

    spin = Spin(M0, T1, T2, Δω / 2π)
    if how == :ideal
        return stfrblochsim(spin, Tfree, Tg, κ * α, κ * β, ϕ, TE)
    elseif how == :grad || how == :rfspoil
        return stfrblochsim_avg(spin, Tfree, Tg, κ * α, κ * β, ϕ, TE, how;
                                kwargs...)
    else
        howerror(how)
    end

end

"""
    stfrblochsim(M0, frac, T1, T2, Δω, τ, κ, Tfree, Tg, α, β, ϕ[, TE];
                 how, kwargs...)

Calculate the STFR signal of a multi-compartment spin using Bloch simulation.

# Arguments
- `M0::Real`: Equilibrium magnetization
- `frac::AbstractArray{<:Real,1}`: Volume fraction of each compartment
- `T1::AbstractArray{<:Real,1}`: Spin-lattice recovery time constants (ms)
- `T2::AbstractArray{<:Real,1}`: Spin-spin recovery time constants (ms)
- `Δω::AbstractArray{<:Real,1}`: Off-resonance frequencies (rad/s)
- `τ::AbstractArray{<:Real,1}`: Residence times (ms); of form
    [τ12, τ13, ..., τ1N, τ21, τ23, ..., τ2N, ...]
- `κ::Real`: Variation in tip angles
- `Tfree::Real`: Free precession time (ms)
- `Tg::Real`: Spoiler gradient time (ms)
- `α::Real`: Tip-down angle (rad)
- `β::Real`: Tip-up angle (rad)
- `ϕ::Real`: Phase of tip-up RF pulse (rad)
- `TE::Real = Tfree / 2`: Echo time (ms)
- `how::Symbol = :ideal`: See docstring for single-compartment `stfrblochsim`

# Return
- `M::Complex{Float64}`: Signal generated from an STFR scan
"""
function stfrblochsim(
    M0::Real,
    frac::AbstractArray{<:Real,1},
    T1::AbstractArray{<:Real,1},
    T2::AbstractArray{<:Real,1},
    Δω::AbstractArray{<:Real,1},
    τ::AbstractArray{<:Real,1},
    κ::Real,
    Tfree::Real,
    Tg::Real,
    α::Real,
    β::Real,
    ϕ::Real,
    TE::Real = Tfree / 2;
    how::Symbol = :ideal,
    kwargs...
)

    spin = SpinMC(M0, frac, T1, T2, Δω / 2π, τ)
    if how == :ideal
        return stfrblochsim(spin, Tfree, Tg, κ * α, κ * β, ϕ, TE)
    elseif how == :grad || how == :rfspoil
        return stfrblochsim_avg(spin, Tfree, Tg, κ * α, κ * β, ϕ, TE, how;
                                kwargs...)
    else
        howerror(how)
    end

end

"""
    stfrblochsim(spin, Tfree, Tg, α, β, ϕ, TE)

Simulate the steady-state signal acquired from STFR, assuming instantaneous
excitations and ideal spoiling.
"""
function stfrblochsim(
    spin::AbstractSpin,
    Tfree::Real,
    Tg::Real,
    α::Real,
    β::Real,
    ϕ::Real,
    TE::Real
)

    (Atd, _) = excitation(spin, 0, α)
    (Atu, _) = excitation(spin, ϕ, -β)
    (Atf, Btf) = freeprecess_nomutate(spin, Tfree)
    (Atg, Btg) = freeprecess_nomutate(spin, Tg)
    S = spoil(spin)
    M = (diagm(0 => trues(size(S, 1))) - Atd * S * Atg * Atu * Atf) \
        (Atd * (S * (Atg * (Atu * Btf))) + Atd * (S * Btg))
    (Ate, Bte) = freeprecess_nomutate(spin, TE)
    M = Ate * M + Bte
    return complex(sum(M[1:3:end]), sum(M[2:3:end]))

end

function stfrblochsim_avg(
    spin::AbstractSpin,
    Tfree::Real,
    Tg::Real,
    α::Real,
    β::Real,
    ϕ::Real,
    TE::Real,
    how::Symbol;
    nspins::Integer = 40,
    ncycles::Real = 1,
    Δθinc::Real = deg2rad(117),
    nTR::Integer = 300
)

    # Pick an arbitrary gradient strength and compute the spatial locations that
    # will provide a uniform ncycles of phase
    gradz = 0.3 # G/cm
    zmax = ncycles * 2π / (GAMMA * gradz * Tg/1000) # cm
    z = (1:nspins)/nspins  * zmax

    # Make nspins copies of the provided spin at different spatial locations
    spins = varyz(spin, z)

    # Compute the STFR signal for each spin
    if how == :grad
        signals = map(spin -> stfrblochsim(spin, Tfree, Tg, α, β, ϕ, TE,
                                           [0,0,gradz]), spins)
    elseif how == :rfspoil
        signals = map(spin -> stfrblochsim(spin, Tfree, Tg, α, β, ϕ, TE,
                                           [0,0,gradz], Δθinc, nTR), spins)
    end

    # Find the average signal
    signal = sum(signals) / nspins

    return signal

end

function varyz(spin::Spin, z)

    map(z -> Spin(spin.M0, spin.T1, spin.T2, spin.Δf, [0,0,z]), z)

end

function varyz(spin::SpinMC, z)

    map(z -> SpinMC(spin.M0, spin.frac, spin.T1, spin.T2, spin.Δf, spin.τ,
                    [0,0,z]), z)

end

"""
    stfrblochsim(spin, Tfree, Tg, α, β, ϕ, TE, grad[, Δθinc, nTR])

Simulate the steady-state signal acquired from STFR, assuming instantaneous
excitations and nonideal spoiling.

The steady-state is a pseudo steady-state if using RF spoiling.
"""
function stfrblochsim(
    spin::AbstractSpin,
    Tfree::Real,
    Tg::Real,
    α::Real,
    β::Real,
    ϕ::Real,
    TE::Real,
    grad::AbstractArray{<:Real}
)

    # Precompute spin dynamics
    Dtd = excitation(spin, 0, α)
    Dte = freeprecess_nomutate(spin, TE)
    Dtr = freeprecess_nomutate(spin, Tfree - TE)
    Dtu = excitation(spin, 0, -β)
    Dtg = freeprecess_nomutate(spin, Tg, grad)

    # Calculate steady-state magnetization immediately following excitation
    (A, B) = combine(Dte, Dtr, Dtu, Dtg, Dtd)
    M = (diagm(0 => trues(size(A, 1))) - A) \ B

    # Calculate steady-state signal at echo time
    M = Dte[1] * M + Dte[2]

    return complex(sum(M[1:3:end]), sum(M[2:3:end]))

end

function stfrblochsim(
    spin::AbstractSpin,
    Tfree::Real,
    Tg::Real,
    α::Real,
    β::Real,
    ϕ::Real,
    TE::Real,
    grad::AbstractArray{<:Real},
    Δθinc::Real,
    nTR::Integer
)

    # Precompute spin dynamics
    Dte = freeprecess_nomutate(spin, TE)
    Dtr = freeprecess_nomutate(spin, Tfree - TE)
    Dtetr = combine(Dte, Dtr)
    Dtg = freeprecess_nomutate(spin, Tg, grad)

    # Initialize RF spoiling parameters
    θ = 0
    Δθ = Δθinc

    # Simulate nTR TR's
    M = spin.M
    for rep = 1:nTR

        (A, B) = excitation(spin, θ, α)
        M = A * M + B
        M = Dtetr[1] * M + Dtetr[2]
        (A, B) = excitation(spin, θ, -β)
        M = A * M + B
        M = Dtg[1] * M + Dtg[2]
        θ += Δθ
        Δθ += Δθinc

    end

    # Calculate signal at echo time
    (A, B) = excitation(spin, θ, α)
    M = A * M + B
    M = Dte[1] * A + Dte[2]

    return complex(sum(M[1:3:end]), sum(M[2:3:end])) * exp(im * θ)

end

# Avoid string interpolation in the case where how is okay
@noinline function howerror(how::Symbol)

    throw(ArgumentError("unsupported how = :$how"))

end

# Versions of freeprecess that doesn't mutate arrays, hence is differentiable
# using Zygote
freeprecess_nomutate(spin::Spin, t::Real) = freeprecess(spin, t)

function freeprecess_nomutate(spin::SpinMC, t::Real)

    E = exp(t * spin.A)
    B = (diagm(0 => trues(size(E, 1))) - E) * spin.Meq
    return (E, B)

end

freeprecess_nomutate(spin::Spin, t::Real, grad::AbstractArray{<:Real,1}) =
    freeprecess(spin, t, grad)

function freeprecess_nomutate(spin::SpinMC, t::Real, grad::AbstractArray{<:Real,1})

    gradfreq = GAMMA * sum(grad .* spin.pos) / 1000 # rad/ms
    ΔA = diagm(1 => repeat([gradfreq, 0, 0], spin.N), # Left-handed rotation
              -1 => repeat([-gradfreq, 0, 0], spin.N))[1:3spin.N,1:3spin.N]
    E = exp(t * (spin.A + ΔA))
    B = (diagm(0 => trues(size(E, 1))) - E) * spin.Meq
    return (E, B)

end
