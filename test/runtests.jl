tests = ["signalmodel", "blochsim"]
for t in tests
    include("$(t).jl")
end
