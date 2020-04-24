function testone1()

    M0 = 1
    T1 = 1000
    T2 = 80
    Δω = 10 * 2π
    κ = 0.9
    Tfree = 8
    Tg = 3
    α = deg2rad(15)
    β = deg2rad(13)
    ϕ = deg2rad(30)
    s1 = stfr(M0, T1, T2, Δω, κ, Tfree, Tg, α, β, ϕ)
    s2 = stfrblochsim(M0, T1, T2, Δω, κ, Tfree, Tg, α, β, ϕ)
    return s1 ≈ s2

end

function testone2()

    M0 = [1.04, 0.8, 1.3]
    T1 = [833, 900, 1100]
    T2 = [76, 90, 100]
    Δω = [10, 0, 30] * 2π
    κ = [0.9, 1, 1.1]
    Tfree = [8, 10]
    Tg = [3, 4]
    α = deg2rad.([15, 14])
    β = deg2rad.([13, 3])
    ϕ = deg2rad.([30, 0])
    s1 = stfr.(M0', T1', T2', Δω', κ', Tfree, Tg, α, β, ϕ)
    s2 = stfrblochsim.(M0', T1', T2', Δω', κ', Tfree, Tg, α, β, ϕ)
    return s1 ≈ s2

end

function testtwo1()

    M0 = 0.8
    ff = 0.15
    T1f = 400
    T1s = 1000
    T2f = 20
    T2s = 80
    Δωf = 17 * 2π
    Δω = 0 * 2π
    κ = 1
    Tfree = 10
    Tg = 5
    α = deg2rad(11)
    β = deg2rad(1)
    ϕ = deg2rad(0)
    frac = [ff, 1-ff]
    T1 = [T1f, T1s]
    T2 = [T2f, T2s]
    Δω2 = [Δωf, 0] .+ Δω
    τ = [Inf, Inf]
    s1 = stfr(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ, Tfree, Tg, α, β, ϕ)
    s2 = stfrblochsim(M0, frac, T1, T2, Δω2, τ, κ, Tfree, Tg, α, β, ϕ)
    return s1 ≈ s2

end

function testtwo2()

    M0 = [1.04, 0.8, 1.3]
    ff = [0.15, 0.10, 0.23]
    T1f = [400, 320, 833]
    T1s = [833, 900, 1100]
    T2f = [20,  16, 24]
    T2s = [76, 90, 100]
    Δωf = [17, 10, 3] * 2π
    Δω = [10, 0, 30] * 2π
    κ = [0.9, 1, 1.1]
    Tfree = [8, 10]
    Tg = [3, 4]
    α = deg2rad.([15, 14])
    β = deg2rad.([13, 3])
    ϕ = deg2rad.([30, 0])
    frac = permutedims([[ff[i], 1 .- ff[i]] for i = 1:length(ff)])
    T1 = permutedims([[T1f[i], T1s[i]] for i = 1:length(T1f)])
    T2 = permutedims([[T2f[i], T2s[i]] for i = 1:length(T2f)])
    Δω2 = permutedims([[Δωf[i], 0] .+ Δω[i] for i = 1:length(Δω)])
    τ = [[Inf, Inf]]
    s1 = stfr.(M0', ff', T1f', T1s', T2f', T2s', Δωf', Δω', κ', Tfree, Tg, α, β, ϕ)
    s2 = stfrblochsim.(M0', frac, T1, T2, Δω2, τ, κ', Tfree, Tg, α, β, ϕ)
    return s1 ≈ s2

end

function testthree1()

    M0 = 1
    ff = 0.15
    fm = 0.1
    T1f = 400
    T1s = 1000
    T1m = 1000
    T2f = 20
    T2s = 80
    T2m = 0.01
    τfs = 100
    τfm = 50
    Δωf = 15 * 2π
    Δω = 0
    κ = 1
    Tfree = 8
    Tg = 3
    α = deg2rad(15)
    β = deg2rad(15)
    ϕ = 0
    TE = 4
    frac = [ff, 1 - ff - fm, fm]
    T1 = [T1f, T1s, T1m]
    T2 = [T2f, T2s, T2m]
    Δω2 = [Δω + Δωf, Δω, Δω]
    τ = [τfs, τfm, τfs * (1 - ff - fm) / ff, Inf, Inf, Inf]
    s1 = stfrblochsim(M0, ff, fm, T1f, T1s, T1m, T2f, T2s, T2m, τfs, τfm, Δωf,
                      Δω, κ, Tfree, Tg, α, β, ϕ, TE)
    s2 = stfrblochsim(M0, frac, T1, T2, Δω2, τ, κ, Tfree, Tg, α, β, ϕ, TE)
    return s1 ≈ s2

end

function testnonideal1()

    M0 = [1.04, 0.8, 1.3]
    T1 = [833, 900, 1100]
    T2 = [76, 90, 100]
    Δω = [10, 0, 30] * 2π
    κ = [0.9, 1, 1.1]
    Tfree = [8, 10]
    Tg = [3, 4]
    α = deg2rad.([15, 14])
    β = deg2rad.([13, 3])
    ϕ = deg2rad.([30, 0])
    s1 = stfrblochsim.(M0', T1', T2', Δω', κ', Tfree, Tg, α, β, ϕ, how = :grad)
    s2 = stfrblochsim.(M0', T1', T2', Δω', κ', Tfree, Tg, α, β, ϕ,
                       how = :rfspoil, Δθinc = 2π, nTR = 2000)
    return s1 ≈ s2

end

function testnonideal2()

    M0 = [1.04, 0.8, 1.3]
    frac = permutedims([[1], [1], [1]])
    T1 = permutedims([[833], [900], [1100]])
    T2 = permutedims([[76], [90], [100]])
    Δω = permutedims([[10], [0], [30]] * 2π)
    τ = permutedims([Int[], Int[], Int[]])
    κ = [0.9, 1, 1.1]
    Tfree = [8, 10]
    Tg = [3, 4]
    α = deg2rad.([15, 14])
    β = deg2rad.([13, 3])
    ϕ = deg2rad.([30, 0])
    s1 = stfrblochsim.(M0', frac, T1, T2, Δω, τ, κ', Tfree, Tg, α, β, ϕ,
                       how = :grad)
    s2 = stfrblochsim.(M0', frac, T1, T2, Δω, τ, κ', Tfree, Tg, α, β, ϕ,
                       how = :rfspoil, Δθinc = 2π, nTR = 2000)
    return s1 ≈ s2

end

function testnoninstant1()

    M0 = [1.04, 0.8, 1.3]
    ff = [0.15, 0.10, 0.23]
    T1f = [400, 320, 833]
    T1s = [833, 900, 1100]
    T2f = [20,  16, 24]
    T2s = [76, 90, 100]
    Δωf = [17, 10, 3] * 2π
    Δω = [10, 0, 30] * 2π
    κ = [0.9, 1, 1.1]
    Tfree = [8, 10]
    Tg = [3, 4]
    α = deg2rad.([15, 14])
    β = deg2rad.([13, 3])
    ϕ = deg2rad.([30, 0])
    frac = permutedims([[ff[i], 1 .- ff[i]] for i = 1:length(ff)])
    T1 = permutedims([[T1f[i], T1s[i]] for i = 1:length(T1f)])
    T2 = permutedims([[T2f[i], T2s[i]] for i = 1:length(T2f)])
    Δω2 = permutedims([[Δωf[i], 0] .+ Δω[i] for i = 1:length(Δω)])
    τ = [[Inf, Inf]]
    s1 = stfr.(M0', ff', T1f', T1s', T2f', T2s', Δωf', Δω', κ', Tfree, Tg, α, β, ϕ)
    s2 = stfrblochsim.(M0', frac, T1, T2, Δω2, τ, κ', Tfree, Tg, α, β, ϕ,
                       rfduration = 0.002, nrf = 11)
    return s1 ≈ s2

end

function testnoninstant2()

    M0 = [1.04, 0.8, 1.3]
    frac = permutedims([[1], [1], [1]])
    T1 = permutedims([[833], [900], [1100]])
    T2 = permutedims([[76], [90], [100]])
    Δω = permutedims([[10], [0], [30]] * 2π)
    τ = permutedims([Int[], Int[], Int[]])
    κ = [0.9, 1, 1.1]
    Tfree = [8, 10]
    Tg = [3, 4]
    α = deg2rad.([15, 14])
    β = deg2rad.([13, 3])
    ϕ = deg2rad.([30, 0])
    s1 = stfrblochsim.(M0', frac, T1, T2, Δω, τ, κ', Tfree, Tg, α, β, ϕ,
                       how = :grad, rfduration = 2, nrf = 11)
    s2 = stfrblochsim.(M0', frac, T1, T2, Δω, τ, κ', Tfree, Tg, α, β, ϕ,
                       how = :rfspoil, rfduration = 2, nrf = 11, Δθinc = 2π,
                       nTR = 2000)
    return s1 ≈ s2

end

function testerror1()

    stfrblochsim(1, 1000, 100, 0, 1, 8, 3, 0, 0, 0, how = :error)
    return false

end

function testerror2()

    stfrblochsim(1, [1], [1000], [100], [0], Int[], 1, 8, 3, 0, 0, 0, how = :error)
    return false

end

@testset "Bloch Simulation" begin

    @testset "Single-Compartment" begin

        @test testone1()
        @test testone2()

    end

    @testset "Two-Compartment" begin

        @test testtwo1()
        @test testtwo2()

    end

    @testset "Three-Compartment" begin

        @test testthree1()

    end

    @testset "Nonideal Spoiling" begin

        @test testnonideal1()
        @test testnonideal2()

    end

    @testset "Non-instantaneous Excitation" begin

        @test testnoninstant1()
        @test testnoninstant2()

    end

    @testset "Errors" begin

        @test_throws ArgumentError testerror1()
        @test_throws ArgumentError testerror2()

    end

end
