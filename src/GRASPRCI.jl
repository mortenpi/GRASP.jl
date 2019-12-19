"""
Functions and types wrapping the `libgrasp-rci` shared library.
"""
module GRASPRCI

function grasp_orbitals_nw()
    @debug "Starting ccall: grasp_orbitals_nw"
    r = ccall( (:grasp_orbitals_nw, libgrasp), Cint, ())
    Int(r)
end

function grasp_orbitals()
    @debug "Starting ccall: grasp_orbitals"
    nw = grasp_orbitals_nw()
    np, nak = Vector{Cint}(undef, nw), Vector{Cint}(undef, nw)
    nw_ref = Ref{Cint}()
    ccall( (:grasp_orbitals, libgrasp),
        Cvoid, (Ptr{Cint}, Ptr{Cint}, Ptr{Cint}),
        nw_ref, np, nak
    )
    return map(args -> RelativisticOrbital(args...), zip(np, nak))
end

end
