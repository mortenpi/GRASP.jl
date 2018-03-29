# Non-relativistic CSF definition lists -- to generate inputs for rcsfgenerate

struct CSFOrbital
    n :: Int
    l :: Int

    function CSFOrbital(n, l)
        @assert n >= 1 && l >= 0
        @assert l < n
        new(n, l)
    end
end

import Base: ==
==(a::CSFOrbital, b::CSFOrbital) = (a.n == b.n) && (a.l == b.l)

Base.isless(a::CSFOrbital, b::CSFOrbital) = (a.n < b.n) || (a.n == b.n && a.l < b.l)

function Base.show(io::IO, co::CSFOrbital)
    write(io, "CSFOrbital(\"")
    print(io, co)
    write(io, "\")")
end

function Base.print(io::IO, co::CSFOrbital)
    write(io, string(co.n), specname(co.l));
    nothing
end


struct CSFDefinition
    orbitals :: Vector{CSFOrbital}
    "Number of allowed excitations."
    nelectrons :: Vector{Int}
    "Number of allowed excitations."
    nexcitations :: Vector{Int}
end
CSFDefinition() = CSFDefinition([], [], [])
function CSFDefinition(csfdef::CSFDefinition)
    CSFDefinition(
        copy(csfdef.orbitals),
        copy(csfdef.nelectrons),
        copy(csfdef.nexcitations)
    )
end

function Base.push!(cdef::CSFDefinition, orbital::CSFOrbital, nelec, nexc)
    # TODO: Sanity checks!
    push!(cdef.orbitals, orbital)
    push!(cdef.nelectrons, nelec)
    push!(cdef.nexcitations, nexc)
    cdef
end

Base.length(cdef::CSFDefinition) = length(cdef.orbitals)

function Base.start(cdef::CSFDefinition)
    # sanity check
    @assert length(cdef.orbitals) == length(cdef.nelectrons) == length(cdef.nexcitations)
    1
end
Base.done(cdef::CSFDefinition, i) = (i > length(cdef.orbitals))
Base.next(cdef::CSFDefinition, i) = (cdef.orbitals[i], cdef.nelectrons[i], cdef.nexcitations[i]), i+1

nelectrons(csf::CSFDefinition) = sum(csf.nelectrons)

function nexcitations(csf_from::CSFDefinition, csf_to::CSFDefinition)
    @assert nelectrons(csf_from) == nelectrons(csf_to)
    let from_orbs = [(orb.n, orb.l) for orb in csf_from.orbitals],
        to_orbs = [(orb.n, orb.l) for orb in csf_to.orbitals]
        @assert length(from_orbs) == length(unique(from_orbs))
        @assert length(to_orbs) == length(unique(to_orbs))
    end

    to_unchecked = fill(true, length(csf_to.orbitals))

    nexcitations = 0
    nremoved = 0 # for sanity checking only
    for (i, orb) in enumerate(csf_from.orbitals)
        # Let's count all the electrons that have been excited into some of
        # the occupied orbitals of csf_from.
        toids = find(o -> o.n == orb.n && o.l == orb.l, csf_to.orbitals)
        @assert length(toids) < 2
        if length(toids) == 0
            nremoved += csf_from.nelectrons[i]
        end
        if length(toids) == 1
            toid = first(toids)
            @assert to_unchecked[toid]
            to_unchecked[toid] = false
            diff = csf_to.nelectrons[toid] - csf_from.nelectrons[i]
            # If diff < 0, electrons have been excited from this orbital
            (diff < 0) && (nremoved += abs(diff)) # for sanity checking
            nexcitations += (diff > 0) ? diff : 0
        end
    end
    for (i, orb) in enumerate(csf_to.orbitals)
        to_unchecked[i] || continue
        nexcitations += csf_to.nelectrons[i]
    end
    @assert nexcitations == nremoved # sanity check
    return nexcitations
end

function Base.print(io::IO, cdef::CSFDefinition)
    for (orb, nelec, nexc) in cdef
        nexc_str = if nexc == nelec
            "*"
        elseif nexc == 0
            "i"
        else
            string(nelec - nexc)
        end
        print(io, orb)
        write(io, "(", string(nelec), ",", nexc_str, ")")
    end
    nothing
end

import Base: *
"""
For two [`CSFDefinition`](@ref) objects, merge the CSF definitions lists together, or thrown
an error if the lists have overlapping orbitals.
"""
function *(a::CSFDefinition, b::CSFDefinition)
    if !isempty(intersect(a.orbitals, b.orbitals))
        isect = join(intersect(a.orbitals, b.orbitals), ",")
        error("""
        CSFDefinitions contain overlapping orbitals ($(isect)).
          1: $(a)
          2: $(b)
        """)
    end

    CSFDefinition(
        vcat(a.orbitals, b.orbitals),
        vcat(a.nelectrons, b.nelectrons),
        vcat(a.nexcitations, b.nexcitations)
    )
end


mutable struct CSFDefinitionList
    cdefs :: Vector{CSFDefinition}
end

function Base.print(io::IO, cdl::CSFDefinitionList)
    first = true
    for cdef in cdl.cdefs
        first ? (first = false) : println(io)
        print(io, cdef)
    end
    nothing
end

"""
    csfdefinition(csf::CSF)

Reduces a `CSF` down to a `CSFDefinition`, by merging the orbitals with same `l` but a
different `j` value.
"""
function csfdefinition(csf::CSF)
    nrelos = Vector{Tuple{Int,Int}}()
    nelecs = Dict{Tuple{Int,Int}, Int}()
    for orb in csf.orbitals
        nl = orb.n, kappa2l(orb.kappa)
        if !(nl in keys(nelecs))
            push!(nrelos, nl)
            nelecs[nl] = 0
        end
        nelecs[nl] += orb.nelectrons
    end

    cd = CSFDefinition()
    for nl in nrelos
        push!(cd, CSFOrbital(nl...), nelecs[nl], 0)
    end
    return cd
end
