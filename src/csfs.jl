const SPECTROSCOPIC_NAMES = "s p d f g h i k l m n o q r t u v" |> split

function specname(l::Integer)
    @assert l >= 0
    SPECTROSCOPIC_NAMES[l+1]
end

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

function Base.show(io::IO, co::CSFOrbital)
    write(io, "CSFOrbital(\"")
    print(io, co)
    write(io, "\")")
end

function Base.print(io::IO, co::CSFOrbital)
    write(io, string(co.n), specname(co.l));
    nothing
end


mutable struct CSFDefinition
    orbitals :: Vector{CSFOrbital}
    "Number of allowed excitations."
    nelectrons :: Vector{Int}
    "Number of allowed excitations."
    nexcitations :: Vector{Int}
    CSFDefinition() = new([],[],[])
end

function Base.push!(cdef::CSFDefinition, orbital::CSFOrbital, nelec, nexc)
    # TODO: Sanity checks!
    push!(cdef.orbitals, orbital)
    push!(cdef.nelectrons, nelec)
    push!(cdef.nexcitations, nexc)
end

Base.length(cdef::CSFDefinition) = length(cdef.orbitals)

function Base.start(cdef::CSFDefinition)
    # sanity check
    @assert length(cdef.orbitals) == length(cdef.nelectrons) == length(cdef.nexcitations)
    1
end
Base.done(cdef::CSFDefinition, i) = (i > length(cdef.orbitals))
Base.next(cdef::CSFDefinition, i) = (cdef.orbitals[i], cdef.nelectrons[i], cdef.nexcitations[i]), i+1

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

    cdef = CSFDefinition()
    cdef.orbitals = vcat(a.orbitals, b.orbitals)
    cdef.nelectrons = vcat(a.nelectrons, b.nelectrons)
    cdef.nexcitations = vcat(a.nexcitations, b.nexcitations)
    return cdef
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
