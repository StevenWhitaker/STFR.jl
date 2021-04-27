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
        new{T1,T2,T3,T4,T5}(Tfree, TE, Tg, rftipdown, rftipup, spoiling)

    end
end

STFRBlochSim(Tfree, Tg, TE, rftipdown, rftipup, spoiling::AbstractSpoiling, nTR::Val = Val(0)) = STFRBlochSim(Tfree, Tg, TE, rftipdown, rftipup, spoiling, nTR, Val(false))
STFRBlochSim(Tfree, Tg, TE, rftipdown, rftipup, nTR::Val, save_transients::Val = Val(false)) = STFRBlochSim(Tfree, Tg, TE, rftipdown, rftipup, IdealSpoiling(), nTR, save_transients)
STFRBlochSim(Tfree, Tg, TE, rftipdown, rftipup) = STFRBlochSim(Tfree, Tg, TE, rftipdown, rftipup, IdealSpoiling(), Val(0), Val(false))
STFRBlochSim(Tfree, Tg, TE, α::Real, rftipup, spoiling, nTR, save_transients) = STFRBlochSim(Tfree, Tg, TE, InstantaneousRF(α), rftipup, spoiling, nTR, save_transients)
STFRBlochSim(Tfree, Tg, TE, rftipdown, β::Real, ϕ::Real, args...) = STFRBlochSim(Tfree, Tg, TE, rftipdown, InstantaneousRF(-β, ϕ), args...)

struct STFRBlochSimWorkspace
    Ate
    Bte
    Atr
    Btr
    Atg
    Btg
    Atd
    Btd
    Atu
    Btu
    BtoM
    bm_workspace::Union{Nothing,<:BlochMcConnellWorkspace} # Maybe? Instead of SpinCollection
end

# Case when nTR = 0
function (scan::STFRBlochSim{<:AbstractRF,<:AbstractRF,<:AbstractSpoiling,0})(
    spin::AbstractSpin,
    workspace::STFRBlochSimWorkspace = STFRBlochSimWorkspace(spin, scan)
)

    rfduration = (duration(scan.rftipdown) + duration(scan.rftipup)) / 2

    excite!(workspace.Atd, workspace.Btd, spin, scan.rftipdown, workspace.ex_workspace)
    freeprecess!(workspace.Atf, workspace.Btf, spin, scan.Tfree - rfduration, workspace.bm_workspace)
    excite!(workspace.Atu, workspace.Btu, spin, scan.rftipup, workspace.ex_workspace)
    freeprecess!(workspace.Atg, workspace.Btg, spin, scan.Tg - rfduration - spoiler_gradient_duration(scan.spoiling), workspace.bm_workspace)
    spoil!(workspace.As, workspace.Bs, spin, scan.spoiling, workspace.bm_workspace)

    combine!(workspace.tmpA1, workspace.tmpB1, workspace.Atf, workspace.Btf, workspace.Atu, workspace.Btu)
    combine!(workspace.tmpA2, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1, workspace.Atg, workspace.Btg)
    combine!(workspace.tmpA1, workspace.tmpB1, workspace.tmpA2, workspace.tmpB2, workspace.Ats, workspace.Bts)
    combine!(workspace.tmpA2, workspace.tmpB2, workspace.tmpA1, workspace.tmpB1, workspace.Atd, workspace.Btd)
    subtract!(workspace.mat, I, workspace.tmpA2)
    copyto!(workspace.vec, workspace.tmpB2)
    F = lu!(workspace.mat)
    ldiv!(F, workspace.vec)
    copyto!(spin.M, workspace.vec)

    freeprecess!(workspace.Atf, workspace.Btf, spin, spin.TE - duration(scan.rftipdown) / 2, workspace.bm_workspace)
    applydynamics!(spin, workspace.tmpB1, workspace.Atf, workspace.Btf)

end
