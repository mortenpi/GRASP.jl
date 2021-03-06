"""
Non-relativistic configurations and CSF definition lists, for generating inputs
to `rcsfgenerate`.
"""
module Configurations

import ..GRASP: angularmomentum, nelectrons, maxelectrons, nexcitations
import ..GRASP: SPECTROSCOPIC_NAMES, CSF, specname
using AtomicLevels: κ2ℓ, κ2j

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

maxelectrons(co::CSFOrbital) = 4*co.l + 2

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

function Base.iterate(cdef::CSFDefinition, s = 1)
    @assert length(cdef.orbitals) == length(cdef.nelectrons) == length(cdef.nexcitations) # sanity check
    s > length(cdef.orbitals) && return nothing
    item = (cdef.orbitals[s], cdef.nelectrons[s], cdef.nexcitations[s])
    return item, s + 1
end

nelectrons(cdef::CSFDefinition) = sum(cdef.nelectrons)

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
        toids = findall(o -> o.n == orb.n && o.l == orb.l, csf_to.orbitals)
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
    for (orb, nelec, orbcoupling, csfcoupling) in csf
        nl = orb.n, κ2ℓ(orb.κ)
        if !(nl in keys(nelecs))
            push!(nrelos, nl)
            nelecs[nl] = 0
        end
        nelecs[nl] += nelec
    end

    cd = CSFDefinition()
    for nl in nrelos
        push!(cd, CSFOrbital(nl...), nelecs[nl], 0)
    end
    return cd
end

const NUMBERS = [Char(48+i) for i=0:9]

function Base.parse(::Type{CSFDefinition}, str)
     # Find all the opening and closing parens first.
    lparens, rparens = Int[], Int[]
    for (i, c) in enumerate(str)
        if c == '('
            (length(lparens) == length(rparens)) || throw(ArgumentError("Unexpected '(' at $i in '$str'"))
            push!(lparens, nextind(str, 0, i))
        elseif c == ')'
            (length(lparens) - length(rparens) == 1) || throw(ArgumentError("Unexpected ')' at $i in '$str'"))
            push!(rparens, nextind(str, 0, i))
        end
    end
    (length(lparens) == length(rparens)) || throw(ArgumentError("Missing final '(' in '$str'"))

    cd = CSFDefinition()
    for i = 1:length(lparens)
        ss_idx = (i == 1) ? 1 : nextind(str, 0, length(str, 1, rparens[i-1]) + 1)
        ss_end = nextind(str, 0, length(str, 1, lparens[i]) - 1)
        ss_nl = str[ss_idx:ss_end]

        ss_idx = nextind(str, 0, length(str, 1, lparens[i]) + 1)
        ss_end = nextind(str, 0, length(str, 1, rparens[i]) - 1)
        ss_el = str[ss_idx:ss_end]

        # parse the n and l value
        ss_nl = strip(ss_nl)

        last_number = 0
        for (i, c) in enumerate(ss_nl)
            c in NUMBERS || break
            last_number = i
        end
        last_number == 0 && throw(ArgumentError("Missing principal value in `$str`"))

        ss_n = ss_nl[1:nextind(ss_nl, 0, last_number)]
        n = parse(Int, ss_n)

        ss_l = ss_nl[nextind(ss_nl, 0, last_number+1):end]
        l = parse_l(ss_l)
        l >= 0 || throw(ArgumentError("Bad spectroscopic orbital `$ss_l` in `$str`"))

        # parse the electron / excitation counts
        ss_el_split = split(ss_el, ',')
        (1 <= length(ss_el_split) <= 2) || throw(ArgumentError("Bad number of terms in ss_el: $ss_el"))
        nelec = parse(Int, ss_el_split[1])
        nexc = if length(ss_el_split) == 1 || ss_el_split[2] == "i"
            0
        elseif ss_el_split[2] == "*"
            nelec
        else
            parse(Int, ss_el_split[2])
        end

        push!(cd, CSFOrbital(n, l), nelec, nexc)
    end

    return cd
end

function parse_l(lstring)
    idx = findfirst(ls -> ls == strip(lstring), SPECTROSCOPIC_NAMES)
    if idx === nothing
        throw(ArgumentError("Unable to parse `$(lstring)`"))
    end
    return idx - 1
end

function parse_j(s)
    if endswith(s,"-")
        @assert parse_l(s[1:end-1]) != 0
        parse_l(s[1:end-1])
    else
        - parse_l(s) - 1
    end
end

end
