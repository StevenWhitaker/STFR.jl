# Test single-compartment STFR for scalar inputs
# Output should be a scalar; size of scalar is ()
function stfr1_size1()

    M0 = 1
    T1 = 1000
    T2 = 80
    Δω = 1 * 2π
    κ = 1.1
    Tfree = 7.55
    Tg = 4.2
    α = 10 * π/180
    β = 10 * π/180
    ϕ = 1 * π/180
    resultsp = stfr(M0, T1, T2, Δω, κ, Tfree, Tg, α, β, ϕ, spoil = true)
    resultun = stfr(M0, T1, T2, Δω, κ, Tfree, Tg, α, β, ϕ, spoil = false)
    return size(resultsp) == () && size(resultun) == ()

end

# Test single-compartment STFR for vector inputs
# Output should be a 2D array
function stfr1_size2()

    M0 = [1,1,1]
    T1 = [1000,1000,1000]
    T2 = [80,80,80]
    Δω = [1,1,1] * 2π
    κ = [1.1,1.1,1.1]
    Tfree = [7.55,7.55]
    Tg = [4.2,4.2]
    α = [10,10] * π/180
    β = [10,10] * π/180
    ϕ = [1,1] * π/180
    resultsp = stfr.(M0', T1', T2', Δω', κ', Tfree, Tg, α, β, ϕ,
                     spoil = true)
    resultun = stfr.(M0', T1', T2', Δω', κ', Tfree, Tg, α, β, ϕ,
                     spoil = false)
    return size(resultsp) == (2,3) && size(resultun) == (2,3)

end

# Test single-compartment STFR for mixed scalar and vector inputs
# Output should be a row vector
function stfr1_size3()

    M0 = [1,1,1]
    T1 = [1000,1000,1000]
    T2 = [80,80,80]
    Δω = [1,1,1] * 2π
    κ = [1.1,1.1,1.1]
    Tfree = 7.55
    Tg = 4.2
    α = 10 * π/180
    β = 10 * π/180
    ϕ = 1 * π/180
    resultsp = stfr.(M0', T1', T2', Δω', κ', Tfree, Tg, α, β, ϕ,
                     spoil = true)
    resultun = stfr.(M0', T1', T2', Δω', κ', Tfree, Tg, α, β, ϕ,
                     spoil = false)
    return size(resultsp) == (1,3) && size(resultun) == (1,3)

end

# Test single-compartment STFR for mixed scalar and vector inputs
# Output should be a vector
function stfr1_size4()

    M0 = 1
    T1 = 1000
    T2 = 80
    Δω = 1 * 2π
    κ = 1.1
    Tfree = [7.55,7.55]
    Tg = [4.2,4.2]
    α = [10,10] * π/180
    β = [10,10] * π/180
    ϕ = [1,1] * π/180
    resultsp = stfr.(M0, T1, T2, Δω, κ, Tfree, Tg, α, β, ϕ, spoil = true)
    resultun = stfr.(M0, T1, T2, Δω, κ, Tfree, Tg, α, β, ϕ, spoil = false)
    return size(resultsp) == (2,) && size(resultun) == (2,)

end

# Test two-compartment STFR for scalar inputs
# Output should be a scalar; size of scalar is ()
function stfr2_size1()

    M0 = 1
    ff = 0.1
    T1f = 400
    T1s = 1000
    T2f = 20
    T2s = 80
    Δωf = 10 * 2π
    Δωs = 1 * 2π
    κ = 1.1
    Tfree = 7.55
    Tg = 4.2
    α = 10 * π/180
    β = 10 * π/180
    ϕ = 1 * π/180
    resultsp = stfr(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δωs, κ, Tfree, Tg, α, β, ϕ,
                    spoil = true)
    resultun = stfr(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δωs, κ, Tfree, Tg, α, β, ϕ,
                    spoil = false)
    return size(resultsp) == () && size(resultun) == ()

end

# Test two-compartment STFR for vector inputs
# Output should be a 2D array
function stfr2_size2()

    M0 = [1,1,1]
    ff = [0.1,0.1,0.1]
    T1f = [400,400,400]
    T1s = [1000,1000,1000]
    T2f = [20,20,20]
    T2s = [80,80,80]
    Δωf = [10,10,10] * 2π
    Δωs = [1,1,1] * 2π
    κ = [1.1,1.1,1.1]
    Tfree = [7.55,7.55]
    Tg = [4.2,4.2]
    α = [10,10] * π/180
    β = [10,10] * π/180
    ϕ = [1,1] * π/180
    resultsp = stfr.(M0', ff', T1f', T1s', T2f', T2s', Δωf', Δωs', κ',
                     Tfree, Tg, α, β, ϕ, spoil = true)
    resultun = stfr.(M0', ff', T1f', T1s', T2f', T2s', Δωf', Δωs', κ',
                     Tfree, Tg, α, β, ϕ, spoil = false)
    return size(resultsp) == (2,3) && size(resultun) == (2,3)

end

# Test two-compartment STFR for mixed scalar and vector inputs
# Output should be a row vector
function stfr2_size3()

    M0 = [1,1,1]
    ff = [0.1,0.1,0.1]
    T1f = [400,400,400]
    T1s = [1000,1000,1000]
    T2f = [20,20,20]
    T2s = [80,80,80]
    Δωf = [10,10,10] * 2π
    Δωs = [1,1,1] * 2π
    κ = [1.1,1.1,1.1]
    Tfree = 7.55
    Tg = 4.2
    α = 10 * π/180
    β = 10 * π/180
    ϕ = 1 * π/180
    resultsp = stfr.(M0', ff', T1f', T1s', T2f', T2s', Δωf', Δωs', κ',
                     Tfree, Tg, α, β, ϕ, spoil = true)
    resultun = stfr.(M0', ff', T1f', T1s', T2f', T2s', Δωf', Δωs', κ',
                     Tfree, Tg, α, β, ϕ, spoil = false)
    return size(resultsp) == (1,3) && size(resultun) == (1,3)

end

# Test two-compartment STFR for mixed scalar and vector inputs
# Output should be a vector
function stfr2_size4()

    M0 = 1
    ff = 0.1
    T1f = 400
    T1s = 1000
    T2f = 20
    T2s = 80
    Δωf = 10 * 2π
    Δωs = 1 * 2π
    κ = 1.1
    Tfree = [7.55,7.55]
    Tg = [4.2,4.2]
    α = [10,10] * π/180
    β = [10,10] * π/180
    ϕ = [1,1] * π/180
    resultsp = stfr.(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δωs, κ, Tfree, Tg, α, β, ϕ,
                     spoil = true)
    resultun = stfr.(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δωs, κ, Tfree, Tg, α, β, ϕ,
                     spoil = false)
    return size(resultsp) == (2,) && size(resultun) == (2,)

end

# Make sure the single-compartment STFR in Julia gives the same results as the
# MATLAB version
# The MATLAB version was validated by reproducing the figures in the papers
function stfr1_mat(matfile::String)

    v = matread(matfile)
    # Check to see if the tissue parameters are vectors; if so, they will be two-
    # dimensional arrays, which are not acceptable; vectorize them
    if length(v["M0"]) > 1
        v["M0"] = v["M0"][:]
        v["T1"] = v["T1"][:]
        v["T2"] = v["T2"][:]
        v["wf"] = v["wf"][:]
        v["kappa"] = v["kappa"][:]
    end
    # Do the same with the scan parameters
    if length(v["Tfree"]) > 1
        v["Tfree"] = v["Tfree"][:]
        v["Tg"] = v["Tg"][:]
        v["alpha"] = v["alpha"][:]
        v["beta"] = v["beta"][:]
        v["phi"] = v["phi"][:]
    end
    result = stfr.(v["M0"]', v["T1"]', v["T2"]', v["wf"]', v["kappa"]',
                   v["Tfree"], v["Tg"], v["alpha"], v["beta"], v["phi"],
                   spoil = v["spoil"])
    # Check if the result is a scalar or not
    if length(result) > 1
        # The MATLAB result will never be a row vector, so if Julia produces a row
        # vector, change it to a normal vector by transposing it
        if length(v["M0"]) > 1 && length(v["Tfree"]) == 1
            result = result'
        end
        result = convert(Array{Complex{Float64},1}, result[:])
        v["result"] = convert(Array{Complex{Float64},1}, v["result"][:])
        err = sum(abs.(result - v["result"]))
    else
        err = abs(result - v["result"])
    end
    return ≈(err, 0.0, atol = 1e-13)

end

# Make sure the two-compartment STFR in Julia gives the same results as the
# MATLAB version
function stfr2_mat(matfile::String)

    v = matread(matfile)
    # Check to see if the tissue parameters are vectors; if so, they will be two-
    # dimensional arrays, which are not acceptable; vectorize them
    if length(v["M0"]) > 1
        v["M0"] = v["M0"][:]
        v["ff"] = v["ff"][:]
        v["T1f"] = v["T1f"][:]
        v["T2f"] = v["T2f"][:]
        v["T1s"] = v["T1s"][:]
        v["T2s"] = v["T2s"][:]
        v["wf"] = v["wf"][:]
        v["wff"] = v["wff"][:]
        v["kappa"] = v["kappa"][:]
    end
    # Do the same with the scan parameters
    if length(v["Tfree"]) > 1
        v["Tfree"] = v["Tfree"][:]
        v["Tg"] = v["Tg"][:]
        v["alpha"] = v["alpha"][:]
        v["beta"] = v["beta"][:]
        v["phi"] = v["phi"][:]
    end
    result = stfr.(v["M0"]', v["ff"]', v["T1f"]', v["T1s"]', v["T2f"]',
                   v["T2s"]', v["wff"]', v["wf"]', v["kappa"]', v["Tfree"],
                   v["Tg"], v["alpha"], v["beta"], v["phi"], spoil = v["spoil"])
    # Check if the result is a scalar or not
    if length(result) > 1
        # The MATLAB result will never be a row vector, so if Julia produces a row
        # vector, change it to a normal vector by transposing it
        if length(v["M0"]) > 1 && length(v["Tfree"]) == 1
            result = result'
        end
        result = convert(Array{Complex{Float64},1}, result[:])
        v["result"] = convert(Array{Complex{Float64},1}, v["result"][:])
        err = sum(abs.(result - v["result"]))
    else
        err = abs(result - v["result"])
    end
    return ≈(err, 0.0, atol = 1e-13)

end

@testset "Signal Model" begin

    @testset "Size" begin

        @test stfr1_size1()
        @test stfr1_size2()
        @test stfr1_size3()
        @test stfr1_size4()
        @test stfr2_size1()
        @test stfr2_size2()
        @test stfr2_size3()
        @test stfr2_size4()

    end

    @testset "Compare to MATLAB" begin

        @test stfr1_mat("compare_matlab/stfr1_test1.mat")
        @test stfr1_mat("compare_matlab/stfr1_test2.mat")
        @test stfr1_mat("compare_matlab/stfr1_test3.mat")
        @test stfr1_mat("compare_matlab/stfr1_test4.mat")
        @test stfr2_mat("compare_matlab/stfr2_test1.mat")
        @test stfr2_mat("compare_matlab/stfr2_test2.mat")
        @test stfr2_mat("compare_matlab/stfr2_test3.mat")
        @test stfr2_mat("compare_matlab/stfr2_test4.mat")

    end

end
