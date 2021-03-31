struct STFRBlochSim
    Tfree
    Tg
    α
    β
    ϕ
    TE
    config # Maybe pass this in separately to function, leaning towards keeping it here though
    # TODO: Can check if TE is valid, i.e., after end of tip-down pulse and before beginning of tip-up pulse
    # Or maybe just throw a warning if TE is during an RF pulse or during Tg
end

struct STFRBlochSimConfig
    rftd # TODO: RF object in BlochSim.jl, includes shape, duration, nrf, etc.
    rftu # TODO: RF object in BlochSim.jl, includes shape, duration, nrf, etc.
    spoiling # TODO: Spoiling object in BlochSim.jl, subtypes include ideal, grad, rf, etc.
    transient # Whether to record transient state, possibly object that contains sampling interval, number of TRs, etc.
end

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

# Different for ideal/gradient spoiling vs rf spoiling
# Different for transient vs not, though can be combined in one function in rf spoiling case

# RF spoiling case
# TODO: Implement
# - is_transient
# - tipdown_duration
# - tipup_duration
# - spoiler_gradient
# - rfspoiling_increment
# - tipdown_pulse
# - tipup_pulse
function (scan::STFRBlochSim)(spin::AbstractSpin, workspace::STFRBlochSimWorkspace)

    transient = is_transient(scan)

    freeprecess!(workspace.Ate, workspace.Bte, spin, scan.TE - tipdown_duration(scan) / 2, workspace.bm_workspace)
    if transient
        freeprecess!(workspace.Atr, workspace.Btr, spin, scan.Tfree - scan.TE - tipup_duration(scan) / 2, workspace.bm_workspace)
    else
        freeprecess!(workspace.Atr, workspace.Btr, spin, scan.Tfree - tipdown_duration(scan) / 2 - tipup_duration(scan) / 2, workspace.bm_workspace)
    end
    freeprecess!(workspace.Atg, workspace.Btg, spin, scan.Tg - tipdown_duration(scan) / 2 - tipup_duration(scan) / 2, spoiler_gradient(scan), workspace.bm_workspace)

    if transient
        Mout = similar(spin.M, 3, numTRs(scan) + 2)
        Mout[:,1] = spin.M
    end
    θ = 0
    Δθ = rfspoiling_increment(scan)
    for rep = 1:numTRs(scan)

        # Tip-down
        excite!(workspace.Atd, workspace.Btd, spin, tipdown_pulse(scan), θ) # TODO: need another workspace?
        applydynamics!(spin, workspace.BtoM, workspace.Atd, workspace.Btd)

        # Free-precession and readout
        if transient
            applydynamics!(spin, workspace.BtoM, workspace.Ate, workspace.Bte)
            Mout[:,rep+1] = spin.M
        end
        applydynamics!(spin, workspace.BtoM, workspace.Atr, workspace.Btr)

        # Tip-up
        excite!(workspace.Atu, workspace.Btu, spin, tipup_pulse(scan), θ + scan.ϕ)
        applydynamics!(spin, workspace.BtoM, workspace.Atu, workspace.Btu)

        # Gradient spoiling
        applydynamics!(spin, workspace.BtoM, workspace.Atg, workspace.Btg)

        # RF spoiling
        θ += Δθ
        Δθ += rfspoiling_increment(scan)

    end

    excite!(workspace.Atd, workspace.Btd, spin, tipdown_pulse(scan), θ)
    applydynamics!(spin, workspace.BtoM, workspace.Atd, workspace.Btd)
    applydynamics!(spin, workspace.BtoM, workspace.Ate, workspace.Bte)
    if transient
        Mout[:,numTRs(scan)+2] = spin.M
    else
        Mout = spin.M
    end

    return Mout

end
