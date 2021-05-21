# STFR.jl

[![action status](https://github.com/StevenWhitaker/STFR.jl/actions/workflows/runtests.yml/badge.svg)](https://github.com/StevenWhitaker/STFR.jl/actions)
[![codecov](https://codecov.io/gh/StevenWhitaker/STFR.jl/branch/main/graph/badge.svg?token=7OQodIVYIV)](https://codecov.io/gh/StevenWhitaker/STFR.jl)

This package provides functionality for simulating the signal acquired from a
small-tip fast recovery (STFR) MRI scan. The function `stfr` implements the
signal equations found in the following papers:

- [J.-F. Nielsen, D. Yoon, and D. C. Noll. Small-tip fast recovery imaging using non-slice-selective tailored tip-up pulses and radiofrequency-spoiling. Mag. Res. Med., 69(3):657-66, March 2013.](https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.24289)
    - For RF spoiled STFR.
- [H. Sun, J. A. Fessler, D. C. Noll, and J.-F. Nielsen. Strategies for improved 3D small-tip fast recovery imaging. Mag. Res. Med., 72(2):389-98, August 2014.](https://onlinelibrary.wiley.com/doi/abs/10.1002/mrm.24947)
    - For RF unspoiled STFR.
    - Note that the paper has an error; the corrected equations were used in
      this package.

The function `stfr` can be used for either a single-compartment model or a
two-compartment model without exchange.

This package also provides a Bloch simulation implementation of STFR
(`STFRBlochSim`) that allows for modeling different RF pulses and spoiling
schemes. (See [BlochSim.jl](https://github.com/StevenWhitaker/BlochSim.jl) for
different subtypes of `AbstractRF` and `AbstractSpoiling`.) In addition,
`STFRBlochSim` can be used for spins with an arbitrary number of compartments
with exchange (see `Spin` and `SpinMC` in
[BlochSim.jl](https://github.com/StevenWhitaker/BlochSim.jl).

## Getting Started

First, type `]` at the Julia REPL to enter the package prompt.
Next, there are two ways to add this package to your Julia environment.
1. **Add via PublicRegistry**: Type
   `registry add https://github.com/StevenWhitaker/PublicRegistry` to add
   PublicRegistry to the list of Julia package registries. Now you can type
   `add STFR` to install the most recent version of STFR as if it were
   registered in Julia's General registry. Furthermore, running `up` to update
   your packages will automatically update STFR as updates become available.
1. **Add via URL**: Type `add https://github.com/StevenWhitaker/STFR.jl#v2.0.0`
   to add STFR v2.0.0 (note that `v2.0.0` can be replaced with whatever version
   is needed). When an update becomes available, you will need to manually add
   the updated version of STFR (e.g., with
   `add https://github.com/StevenWhitaker/STFR.jl#v2.0.1`).
After the package has been installed, hit backspace to return to the normal
Julia prompt, and then type `using STFR` to load the package.

## Examples

Below are some concrete examples of how to use this package.

```julia
julia> using STFR

julia> sig = stfr(1, 1000, 100, 3.75 * 2π, 1, 8, 3, π/12, π/12, 0, 4) # Spoiled STFR
0.1587855247046399 - 0.01500965127299392im

julia> stfr(1, 1000, 100, 3.75 * 2π, 1, 8, 3, π/12, π/12, 0, 4, Val(false)) # No RF spoiling
0.009194809492086849 + 0.16021302758102335im

julia> spin = Spin(1, 1000, 100, 3.75)
Spin{Float64}:
 M = Magnetization(0.0, 0.0, 1.0)
 M0 = 1.0
 T1 = 1000.0 ms
 T2 = 100.0 ms
 Δf = 3.75 Hz
 pos = Position(0.0, 0.0, 0.0) cm

julia> stfr! = STFRBlochSim(8, 3, 4, π/12, π/12, 0) # Create an object to simulate an STFR scan
Small-Tip Fast Recovery (STFR) Bloch Simulation:
 Tfree = 8 ms
 Tg = 3 ms
 TE = 4 ms
 rftipdown (tip-down excitation pulse) = Instantaneous RF pulse with eltype Float64:
 α = 0.2617993877991494 rad
 θ = 0.0 rad
 rftipup (tip-up RF pulse) = Instantaneous RF pulse with eltype Float64:
 α = -0.2617993877991494 rad
 θ = 0.0 rad
 spoiling = IdealSpoiling()
 steady-state

julia> stfr!(spin) # Simulate a steady-state STFR scan applied to the given spin

julia> spin.M # Steady-state magnetization
Magnetization vector with eltype Float64:
 Mx = 0.15878552470463958
 My = -0.015009651272993887
 Mz = 0.6210482688962908

julia> signal(spin) ≈ sig
true
```
