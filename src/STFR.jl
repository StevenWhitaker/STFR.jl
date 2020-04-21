module STFR

using BlochSim
using LinearAlgebra: Diagonal

include("signalmodel.jl")
include("blochsim.jl")

export stfr
export stfrblochsim

end
