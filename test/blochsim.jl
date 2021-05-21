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

    spin = Spin(M0, T1, T2, Δω / 2π)
    stfrblochsim! = STFRBlochSim(Tfree, Tg, Tfree / 2, κ*α, κ*β, ϕ)
    stfrblochsim!(spin)

    return s1 ≈ signal(spin)

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

    s2 = Matrix{ComplexF64}(undef, length(Tfree), length(M0))
    workspace = STFRBlochSimWorkspace(Spin{Float64}, STFRBlochSim{<:Real,InstantaneousRF,InstantaneousRF,IdealSpoiling})
    for j = 1:length(M0), i = 1:length(Tfree)
        spin = Spin(M0[j], T1[j], T2[j], Δω[j] / 2π)
        stfrblochsim! = STFRBlochSim(Tfree[i], Tg[i], Tfree[i] / 2, κ[j] * α[i], κ[j] * β[i], ϕ[i])
        stfrblochsim!(spin, workspace)
        s2[i,j] = signal(spin)
    end

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
    s1 = stfr(M0, ff, T1f, T1s, T2f, T2s, Δωf, Δω, κ, Tfree, Tg, α, β, ϕ)

    frac = (ff, 1-ff)
    T1 = (T1f, T1s)
    T2 = (T2f, T2s)
    Δω2 = (Δωf + Δω, Δω)
    τ = (Inf, Inf)
    spin = SpinMC(M0, frac, T1, T2, Δω2 ./ 2π, τ)
    stfrblochsim! = STFRBlochSim(Tfree, Tg, Tfree / 2, κ*α, κ*β, ϕ)
    stfrblochsim!(spin)

    return s1 ≈ signal(spin)

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

    s2 = Matrix{ComplexF64}(undef, length(Tfree), length(M0))
    workspace = STFRBlochSimWorkspace(SpinMC{Float64,2}, STFRBlochSim{<:Real,InstantaneousRF,InstantaneousRF,IdealSpoiling})
    for j = 1:length(M0), i = 1:length(Tfree)
        spin = SpinMC(M0[j], (ff[j], 1 - ff[j]), (T1f[j], T1s[j]), (T2f[j], T2s[j]), (Δωf[j] + Δω[j], Δω[j]) ./ 2π, (Inf, Inf))
        stfrblochsim! = STFRBlochSim(Tfree[i], Tg[i], Tfree[i] / 2, κ[j] * α[i], κ[j] * β[i], ϕ[i])
        stfrblochsim!(spin, workspace)
        s2[i,j] = signal(spin)
    end

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
    s1 = Matrix{ComplexF64}(undef, length(Tfree), length(M0))
    s2 = similar(s1)
    s3 = similar(s1)
    w1 = STFRBlochSimWorkspace(Spin{Float64}, STFRBlochSim{<:Real,InstantaneousRF,InstantaneousRF,GradientSpoiling})
    w2 = STFRBlochSimWorkspace(Spin{Float64}, STFRBlochSim{<:Real,InstantaneousRF,InstantaneousRF,RFSpoiling})
    w3 = STFRBlochSimWorkspace(Spin{Float64}, STFRBlochSim{<:Real,InstantaneousRF,InstantaneousRF,RFandGradientSpoiling})
    spoil1 = GradientSpoiling(0, 0, 0, 2)
    spoil2 = RFSpoiling(2π)
    spoil3 = RFandGradientSpoiling(spoil1, spoil2)
    nTR = Val(2000)
    for j = 1:length(M0), i = 1:length(Tfree)
        spin = Spin(M0[j], T1[j], T2[j], Δω[j] / 2π, Position(1, 1, 1))
        stfrblochsim! = STFRBlochSim(Tfree[i], Tg[i], Tfree[i] / 2, κ[j] * α[i], κ[j] * β[i], ϕ[i], spoil1)
        stfrblochsim!(spin, w1)
        s1[i,j] = signal(spin)
        spin.M.x = 0; spin.M.y = 0; spin.M.z = M0[j]
        stfrblochsim! = STFRBlochSim(Tfree[i], Tg[i], Tfree[i] / 2, κ[j] * α[i], κ[j] * β[i], ϕ[i], spoil2, nTR)
        stfrblochsim!(spin, w2)
        s2[i,j] = signal(spin)
        spin.M.x = 0; spin.M.y = 0; spin.M.z = M0[j]
        stfrblochsim! = STFRBlochSim(Tfree[i], Tg[i], Tfree[i] / 2, κ[j] * α[i], κ[j] * β[i], ϕ[i], spoil3, nTR)
        stfrblochsim!(spin, w3)
        s3[i,j] = signal(spin)
    end

    return s1 ≈ s2 ≈ s3

end

function testnonideal2()

    M0 = [1.04, 0.8, 1.3]
    frac = [(0.1, 0.9), (0.15, 0.85), (0.2, 0.8)]
    T1 = [(300, 833), (400, 900), (800, 1100)]
    T2 = [(10, 76), (15, 90), (20, 100)]
    Δf = [(0, 0), (10, 0), (30, 10)]
    τ = [(20, 100), (30, 130), (70, 70)]
    κ = [0.9, 1, 1.1]
    Tfree = [8, 10]
    Tg = [3, 4]
    α = deg2rad.([15, 14])
    β = deg2rad.([13, 3])
    ϕ = deg2rad.([30, 0])
    s1 = Matrix{ComplexF64}(undef, length(Tfree), length(M0))
    s2 = similar(s1)
    s3 = similar(s1)
    w1 = STFRBlochSimWorkspace(SpinMC{Float64,2}, STFRBlochSim{<:Real,InstantaneousRF,InstantaneousRF,GradientSpoiling})
    w2 = STFRBlochSimWorkspace(SpinMC{Float64,2}, STFRBlochSim{<:Real,InstantaneousRF,InstantaneousRF,RFSpoiling})
    w3 = STFRBlochSimWorkspace(SpinMC{Float64,2}, STFRBlochSim{<:Real,InstantaneousRF,InstantaneousRF,RFandGradientSpoiling})
    spoil1 = GradientSpoiling(0, 0, 0, 2)
    spoil2 = RFSpoiling(2π)
    spoil3 = RFandGradientSpoiling(spoil1, spoil2)
    nTR = Val(2000)
    for j = 1:length(M0), i = 1:length(Tfree)
        spin = SpinMC(M0[j], frac[j], T1[j], T2[j], Δf[j], τ[j], Position(1, 1, 1))
        stfrblochsim! = STFRBlochSim(Tfree[i], Tg[i], Tfree[i] / 2, κ[j] * α[i], κ[j] * β[i], ϕ[i], spoil1)
        stfrblochsim!(spin, w1)
        s1[i,j] = signal(spin)
        for k = 1:2 spin.M[k].x = 0; spin.M[k].y = 0; spin.M[k].z = frac[j][k] * M0[j] end
        stfrblochsim! = STFRBlochSim(Tfree[i], Tg[i], Tfree[i] / 2, κ[j] * α[i], κ[j] * β[i], ϕ[i], spoil2, nTR)
        stfrblochsim!(spin, w2)
        s2[i,j] = signal(spin)
        for k = 1:2 spin.M[k].x = 0; spin.M[k].y = 0; spin.M[k].z = frac[j][k] * M0[j] end
        stfrblochsim! = STFRBlochSim(Tfree[i], Tg[i], Tfree[i] / 2, κ[j] * α[i], κ[j] * β[i], ϕ[i], spoil3, nTR)
        stfrblochsim!(spin, w3)
        s3[i,j] = signal(spin)
    end

    return s1 ≈ s2 ≈ s3

end

function testnonideal3()

    Tfree = 10
    Tg = 3
    TE = Tfree / 2
    α = deg2rad(15)
    β = deg2rad(13)
    ϕ = 0
    (M0, T1, T2, Δf) = (1, 1000, 100, 0)
    spin = Spin(M0, T1, T2, Δf)
    stfrblochsim! = STFRBlochSim(Tfree, Tg, TE, α, β, ϕ)
    stfrblochsim!(spin)
    M1 = copy(spin.M)

    grad = Gradient(0, 0, 0.3)
    zmax = 2π / (GAMMA * grad.z * Tg / 1000)
    nz = 100
    z = (1:nz)/nz * zmax
    spoil = RFandGradientSpoiling(grad, Tg)
    stfrblochsim! = STFRBlochSim(Tfree, Tg, TE, α, β, ϕ, spoil, Val(2000))
    workspace = STFRBlochSimWorkspace(spin, stfrblochsim!)
    M2 = Magnetization(0.0, 0.0, 0.0)
    foreach(z) do z
        spin.M.x = 0; spin.M.y = 0; spin.M.z = M0; spin.pos.z = z
        stfrblochsim!(spin, workspace)
        add!(M2, spin.M)
    end
    div!(M2, nz)

    return isapprox(M1, M2; atol = 1e-4)

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
    frac = [(ff[i], 1 - ff[i]) for i = 1:length(ff)]
    T1 = collect(zip(T1f, T1s))
    T2 = collect(zip(T2f, T2s))
    Δf = [(Δωf[i] + Δω[i], Δω[i]) ./ 2π for i = 1:length(Δωf)]
    τ = (Inf, Inf)
    Δt = 0.0002
    s1 = stfr.(M0', ff', T1f', T1s', T2f', T2s', Δωf', Δω', κ', Tfree, Tg, α, β, ϕ)

    s2 = Matrix{ComplexF64}(undef, length(Tfree), length(M0))
    workspace = STFRBlochSimWorkspace(SpinMC{Float64,2}, STFRBlochSim{<:Real,RF,RF,IdealSpoiling})
    for j = 1:length(M0), i = 1:length(Tfree)
        spin = SpinMC(M0[j], frac[j], T1[j], T2[j], Δf[j], τ)
        rftipdown = RF(fill(κ[j] * α[i] * 1000 / Δt / GAMMA / 3, 3), Δt)
        rftipup = RF(fill(κ[j] * -β[i] * 1000 / Δt / GAMMA / 3, 3), Δt, ϕ[i])
        stfrblochsim! = STFRBlochSim(Tfree[i], Tg[i], Tfree[i] / 2, rftipdown, rftipup)
        stfrblochsim!(spin, workspace)
        s2[i,j] = signal(spin)
    end

    return isapprox(s1, s2, atol = 1e-5)

end

function testnoninstant2()

    M0 = [1.04, 0.8, 1.3]
    frac = [(0.1, 0.9), (0.15, 0.85), (0.2, 0.8)]
    T1 = [(300, 833), (400, 900), (800, 1100)]
    T2 = [(10, 76), (15, 90), (20, 100)]
    Δf = [(0, 0), (10, 0), (30, 10)]
    τ = [(20, 100), (30, 130), (70, 70)]
    κ = [0.9, 1, 1.1]
    Tfree = [8, 10]
    Tg = [3, 4]
    α = deg2rad.([15, 14])
    β = deg2rad.([13, 3])
    ϕ = deg2rad.([30, 0])
    Δt = 0.5
    s1 = Matrix{ComplexF64}(undef, length(Tfree), length(M0))
    s2 = similar(s1)
    w1 = STFRBlochSimWorkspace(SpinMC{Float64,2}, STFRBlochSim{<:Real,RF,RF,GradientSpoiling})
    w2 = STFRBlochSimWorkspace(SpinMC{Float64,2}, STFRBlochSim{<:Real,RF,RF,RFandGradientSpoiling})
    spoil1 = GradientSpoiling(1, 2, 3, 1)
    spoil2 = RFandGradientSpoiling(spoil1, 2π)
    nTR = Val(2000)
    for j = 1:length(M0), i = 1:length(Tfree)
        spin = SpinMC(M0[j], frac[j], T1[j], T2[j], Δf[j], τ[j], Position(1, 1, 1))
        rftipdown = RF(fill(κ[j] * α[i] * 1000 / Δt / GAMMA / 3, 3), Δt)
        rftipup = RF(fill(κ[j] * -β[i] * 1000 / Δt / GAMMA / 3, 3), Δt, ϕ[i])
        stfrblochsim! = STFRBlochSim(Tfree[i], Tg[i], Tfree[i] / 2, rftipdown, rftipup, spoil1)
        stfrblochsim!(spin, w1)
        s1[i,j] = signal(spin)
        for k = 1:2 spin.M[k].x = 0; spin.M[k].y = 0; spin.M[k].z = frac[j][k] * M0[j] end
        stfrblochsim! = STFRBlochSim(Tfree[i], Tg[i], Tfree[i] / 2, rftipdown, rftipup, spoil2, nTR)
        stfrblochsim!(spin, w2)
        s2[i,j] = signal(spin)
    end

    return s1 ≈ s2

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

    @testset "Nonideal Spoiling" begin

        @test testnonideal1()
        @test testnonideal2()
        @test testnonideal3()

    end

    @testset "Non-instantaneous Excitation" begin

        @test testnoninstant1()
        @test testnoninstant2()

    end

end
