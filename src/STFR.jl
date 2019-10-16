module STFR

using LinearAlgebra: Diagonal
using BlochSim

include("signalmodel.jl")
include("blochsim.jl")

export stfr
export stfrblochsim

end
