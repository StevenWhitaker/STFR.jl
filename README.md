# STFR.jl

[![Build Status](https://travis-ci.org/StevenWhitaker/STFR.jl.svg?branch=master)](https://travis-ci.org/StevenWhitaker/STFR.jl)
[![codecov](https://codecov.io/gh/StevenWhitaker/STFR.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/StevenWhitaker/STFR.jl)

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
two-compartment model without exchange. The function `stfrblochsim` computes
the STFR signal using Bloch simulation and can be used for an arbitrary number
of compartments with exchange.

## Getting Started
At the Julia REPL, type `]` to enter the package prompt. Then type
`add https://github.com/StevenWhitaker/STFR.jl#v1.0.0` to add STFR v1.0.0
(note that `v1.0.0` can be replaced with whatever version is needed). Hit
backspace to return to the normal Julia prompt, and then type `using STFR` to
load the package.
