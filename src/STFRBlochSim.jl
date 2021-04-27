struct STFRBlochSim{T1<:AbstractRF,T2<:AbstractRF,T3<:AbstractSpoiling,nTR,save_transients}
    Tfree::Float64
    Tg::Float64
    TE::Float64
    rftipdown::T1
    rftipup::T2
    spoiling::T3

    function STFRBlochSim(
        Tfree::Real,
        Tg::Real,
        TE::Real,
        rftipdown::T1,
        rftipup::T2,
        spoiling::T3,
        nTR::Val{T4},
        save_transients::Val{T5}
    ) where {T1<:AbstractRF,T2<:AbstractRF,T3<:AbstractSpoiling,T4,T5}

        Tfree >= TE + duration(rftipup) / 2 ||
            error("Tfree must be greater than or equal to TE + duration(rftipup) / 2")
        TE >= duration(rftipdown) / 2 || error("TE must not be during the tip-down pulse")
        Tg >= duration(rftipup) / 2 + spoiler_gradient_duration(spoiling) + duration(rftipdown) / 2 ||
            error("Tg must be greater than or equal to duration(rftipup) / 2 + " *
                  "duration(rftipdown) / 2 + spoiler_gradient_duration")
        (T4 isa Int && T4 >= 0) || error("nTR must be a nonnegative Int")
        T3 <: Union{<:RFSpoiling,<:RFandGradientSpoiling} && (T4 > 0 ||
            error("nTR must be positive when simulating RF spoiling"))
        T5 isa Bool || error("save_transients must be a Bool")
        T4 == 0 && T5 &&
            @warn("save_transients is true, but nTR = 0; no transients will be saved")
        new{T1,T2,T3,T4,T5}(Tfree, Tg, TE, rftipdown, rftipup, spoiling)

    end
end

STFRBlochSim(Tfree, Tg, TE, rftipdown, rftipup, spoiling::AbstractSpoiling, nTR::Val = Val(0)) = STFRBlochSim(Tfree, Tg, TE, rftipdown, rftipup, spoiling, nTR, Val(false))
STFRBlochSim(Tfree, Tg, TE, rftipdown, rftipup, nTR::Val, save_transients::Val = Val(false)) = STFRBlochSim(Tfree, Tg, TE, rftipdown, rftipup, IdealSpoiling(), nTR, save_transients)
STFRBlochSim(Tfree, Tg, TE, rftipdown, rftipup) = STFRBlochSim(Tfree, Tg, TE, rftipdown, rftipup, IdealSpoiling(), Val(0), Val(false))
STFRBlochSim(Tfree, Tg, TE, α::Real, rftipup, spoiling::AbstractSpoiling, nTR, save_transients) = STFRBlochSim(Tfree, Tg, TE, InstantaneousRF(α), rftipup, spoiling, nTR, save_transients)
STFRBlochSim(Tfree, Tg, TE, rftipdown, β::Real, ϕ::Real, args...) = STFRBlochSim(Tfree, Tg, TE, rftipdown, InstantaneousRF(-β, ϕ), args...)

struct STFRBlochSimWorkspace{T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,T12,T13,T14}
    Atd::T1
    Btd::T2
    Ate::T3
    Bte::T4
    Atf::T3
    Btf::T4
    Atu::T5
    Btu::T6
    Atg::T3
    Btg::T4
    As::T7
    Bs::T8
    tmpA1::T9
    tmpB1::T4
    tmpA2::T9
    tmpB2::T4
    mat::T10
    vec::T11
    bm_workspace::T12
    td_workspace::T13
    tu_workspace::T14
end

function STFRBlochSimWorkspace(
    spin::AbstractSpin,
    scan::STFRBlochSim,
    bm_workspace = spin isa Spin ? nothing : BlochMcConnellWorkspace(spin)
)

    STFRBlochSimWorkspace(typeof(spin), typeof(scan), bm_workspace)

end

function STFRBlochSimWorkspace(
    spin::Union{Type{Spin{T}},Type{SpinMC{T,N}}},
    scan::Type{<:STFRBlochSim{T1,T2,T3}},
    bm_workspace = spin <: Spin ? nothing : BlochMcConnellWorkspace(spin)
) where {T,N,T1,T2,T3}

    if T1 <: InstantaneousRF
        Atd = ExcitationMatrix{T}()
        Btd = nothing
        td_workspace = nothing
    else
        if spin <: Spin
            Atd = BlochMatrix{T}()
            Btd = Magnetization{T}()
        else
            Atd = BlochMcConnellMatrix{T}(N)
            Btd = MagnetizationMC{T}(N)
        end
        td_workspace = ExcitationWorkspace(spin, bm_workspace)
    end
    if T2 <: InstantaneousRF
        Atu = ExcitationMatrix{T}()
        Btu = nothing
        tu_workspace = nothing
    else
        if spin <: Spin
            Atu = BlochMatrix{T}()
            Btu = Magnetization{T}()
        else
            Atu = BlochMcConnellMatrix{T}(N)
            Btu = MagnetizationMC{T}(N)
        end
        tu_workspace = ExcitationWorkspace(spin, bm_workspace)
    end
    if T3 <: IdealSpoiling
        As = idealspoiling
        Bs = nothing
    elseif T3 <: RFSpoiling
        As = nothing
        Bs = nothing
    elseif spin <: Spin
        As = FreePrecessionMatrix{T}()
        Bs = Magnetization{T}()
    else
        As = BlochMcConnellMatrix{T}(N)
        Bs = MagnetizationMC{T}(N)
    end
    if spin <: Spin
        Ate = FreePrecessionMatrix{T}()
        Bte = Magnetization{T}()
        Atf = FreePrecessionMatrix{T}()
        Btf = Magnetization{T}()
        Atg = FreePrecessionMatrix{T}()
        Btg = Magnetization{T}()
        tmpA1 = BlochMatrix{T}()
        tmpB1 = Magnetization{T}()
        tmpA2 = BlochMatrix{T}()
        tmpB2 = Magnetization{T}()
        mat = Matrix{T}(undef, 3, 3)
        vec = Vector{T}(undef, 3)
    else
        Ate = BlochMcConnellMatrix{T}(N)
        Bte = MagnetizationMC{T}(N)
        Atf = BlochMcConnellMatrix{T}(N)
        Btf = MagnetizationMC{T}(N)
        Atg = BlochMcConnellMatrix{T}(N)
        Btg = MagnetizationMC{T}(N)
        tmpA1 = BlochMcConnellMatrix{T}(N)
        tmpB1 = MagnetizationMC{T}(N)
        tmpA2 = BlochMcConnellMatrix{T}(N)
        tmpB2 = MagnetizationMC{T}(N)
        mat = Matrix{T}(undef, 6, 6)
        vec = Vector{T}(undef, 6)
    end
    STFRBlochSimWorkspace(Atd, Btd, Ate, Bte, Atf, Btf, Atu, Btu, Atg, Btg, As,
                          Bs, tmpA1, tmpB1, tmpA2, tmpB2, mat, vec,
                          bm_workspace, td_workspace, tu_workspace)

end

# Case when nTR = 0
function (scan::STFRBlochSim{<:AbstractRF,<:AbstractRF,<:AbstractSpoiling,0})(
    spin::AbstractSpin,
    workspace::STFRBlochSimWorkspace = STFRBlochSimWorkspace(spin, scan)
)

    rfduration = (duration(scan.rftipdown) + duration(scan.rftipup)) / 2

    excite!(workspace.Atd, workspace.Btd, spin, scan.rftipdown, workspace.td_workspace)
    freeprecess!(workspace.Atf, workspace.Btf, spin, scan.Tfree - rfduration, workspace.bm_workspace)
    excite!(workspace.Atu, workspace.Btu, spin, scan.rftipup, workspace.tu_workspace)
    freeprecess!(workspace.Atg, workspace.Btg, spin, scan.Tg - rfduration - spoiler_gradient_duration(scan.spoiling), workspace.bm_workspace)
    spoil!(workspace.As, workspace.Bs, spin, scan.spoiling, workspace.bm_workspace)

    combine!(workspace.tmpA1, workspace.tmpB1, workspace.Atf, workspace.Btf, workspace.Atu, workspace.Btu)
    combine!(workspace.tmpA2, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1, workspace.Atg, workspace.Btg)
    combine!(workspace.tmpA1, workspace.tmpB1, workspace.tmpA2, workspace.tmpB2, workspace.As, workspace.Bs)
    combine!(workspace.tmpA2, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1, workspace.Atd, workspace.Btd)
    subtract!(workspace.mat, I, workspace.tmpA2)
    copyto!(workspace.vec, workspace.tmpB2)
    F = lu!(workspace.mat)
    ldiv!(F, workspace.vec)
    copyto!(spin.M, workspace.vec)

    freeprecess!(workspace.Ate, workspace.Bte, spin, scan.TE - duration(scan.rftipdown) / 2, workspace.bm_workspace)
    applydynamics!(spin, workspace.tmpB1, workspace.Ate, workspace.Bte)

end

# Case when nTR > 0
function (scan::STFRBlochSim{<:AbstractRF,<:AbstractRF,T,nTR,save})(
    spin::AbstractSpin,
    workspace::STFRBlochSimWorkspace = STFRBlochSimWorkspace(spin, scan)
) where {T,nTR,save}

    rftd = scan.rftipdown
    rftu = scan.rftipup
    rfspoiling = T <: Union{<:RFSpoiling,<:RFandGradientSpoiling}
    rfduration = (duration(rftd) + duration(rftu)) / 2

    if rfspoiling
        rftd isa RF && (rftd.Δθ[] = rftd.Δθ_initial)
        rftu isa RF && (rftu.Δθ[] = rftu.Δθ_initial)
        Δθinc = rfspoiling_increment(scan.spoiling)
        θ = zero(Δθinc) # For knowing how much phase to remove when recording signal
        Δθ = Δθinc
    else
        excite!(workspace.Atd, workspace.Btd, spin, rftd, workspace.td_workspace)
        excite!(workspace.Atu, workspace.Btu, spin, rftu, workspace.tu_workspace)
    end

    if save
        M = Vector{typeof(spin.M)}(undef, nTR)
        freeprecess!(workspace.Atf, workspace.Btf, spin, scan.Tfree - scan.TE - duration(rftu) / 2, workspace.bm_workspace)
    else
        freeprecess!(workspace.Atf, workspace.Btf, spin, scan.Tfree - rfduration, workspace.bm_workspace)
    end
    freeprecess!(workspace.Atg, workspace.Btg, spin, scan.Tg - rfduration - spoiler_gradient_duration(scan.spoiling), workspace.bm_workspace)
    spoil!(workspace.As, workspace.Bs, spin, scan.spoiling, workspace.bm_workspace)
    combine!(workspace.tmpA1, workspace.tmpB1, workspace.Atg, workspace.Btg, workspace.As, workspace.Bs)
    freeprecess!(workspace.Ate, workspace.Bte, spin, scan.TE - duration(rftd) / 2, workspace.bm_workspace)

    for rep = 1:nTR-1

        rfspoiling && excite!(workspace.Atd, workspace.Btd, spin, rftd, workspace.td_workspace)
        applydynamics!(spin, workspace.tmpB2, workspace.Atd, workspace.Btd)
        if save
            applydynamics!(spin, workspace.tmpB2, workspace.Ate, workspace.Bte)
            M[rep] = copy(spin.M)
            if rfspoiling
                modulation = exp(im * θ)
                if spin isa Spin
                    tmp = signal(spin.M) * modulation
                    M[rep].x = real(tmp)
                    M[rep].y = imag(tmp)
                else
                    for i = 1:spin.N
                        tmp = signal(spin.M[i]) * modulation
                        M[rep][i].x = real(tmp)
                        M[rep][i].y = imag(tmp)
                    end
                end
            end
        end
        applydynamics!(spin, workspace.tmpB2, workspace.Atf, workspace.Btf)
        rfspoiling && excite!(workspace.Atu, workspace.Btu, spin, rftu, workspace.tu_workspace)
        applydynamics!(spin, workspace.tmpB2, workspace.Atu, workspace.Btu)
        applydynamics!(spin, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1)

        if rfspoiling
            if rftd isa InstantaneousRF
                rftd = InstantaneousRF(rftd.α, rftd.θ + Δθ)
            else
                rftd.Δθ[] += Δθ
            end
            if rftu isa InstantaneousRF
                rftu = InstantaneousRF(rftu.α, rftu.θ + Δθ)
            else
                rftu.Δθ[] += Δθ
            end
            θ += Δθ
            Δθ += Δθinc
        end

    end

    rfspoiling && excite!(workspace.Atd, workspace.Btd, spin, rftd, workspace.td_workspace)
    applydynamics!(spin, workspace.tmpB2, workspace.Atd, workspace.Btd)
    applydynamics!(spin, workspace.tmpB2, workspace.Ate, workspace.Bte)
    if rfspoiling
        modulation = exp(im * θ)
        if spin isa Spin
            tmp = signal(spin.M) * modulation
            spin.M.x = real(tmp)
            spin.M.y = imag(tmp)
        else
            for i = 1:spin.N
                tmp = signal(spin.M[i]) * modulation
                spin.M[i].x = real(tmp)
                spin.M[i].y = imag(tmp)
            end
        end
    end
    if save
        M[nTR] = copy(spin.M)
        return M
    end

end
