module STFR

using LinearAlgebra: diagm
using BlochSim

include("signalmodel.jl")
include("blochsim.jl")

export stfr
export stfrblochsim

end
