module CSLGeneration
using AtomicLevels
using HalfIntegers

export write_rcsfinput

struct CSFDefConfiguration
    c :: Configuration{Orbital{Int}}
    minelecs :: Vector{Int}

    function CSFDefConfiguration(c::Configuration, minelecs)
        length(c) == length(minelecs) || throw(ArgumentError("Configuration and min. electron needs to be of the same length"))
        for i = 1:length(c)
            o, nelec, _ = c[i]
            nelec >= minelecs[i] || throw(ArgumentError("Invalid min.elec ($(minelecs[i])) for $o. Must be <= $(nelec)"))
        end
        new(c, minelecs)
    end
end

CSFDefConfiguration(c::Configuration) = CSFDefConfiguration(c, [num_electrons(c, o) for (o, _, _) in c])

Base.length(cdefconf::CSFDefConfiguration) = length(cdefconf.c)

function Base.iterate(cdefconf::CSFDefConfiguration, i = 1)
    if i <= length(cdefconf)
        ((cdefconf.c[i][1], cdefconf.c[i][2], cdefconf.minelecs[i]), i+1)
    else
        nothing
    end
end

struct CSFDefinition
    configurations :: Vector{CSFDefConfiguration}
    maxorbitals :: Vector{Orbital{Int}}
    nsubstitutions :: Int
    jvalues :: Union{Tuple{HalfInteger,HalfInteger},Nothing}
end

function CSFDefinition(configurations::Vector{CSFDefConfiguration}, maxorbitals::Vector{Orbital{Int}}, nsubstitutions::Int)
    CSFDefinition(configurations, maxorbitals, nsubstitutions, nothing)
end

function largestcore(cs)
    cores = AtomicLevels.get_noble_core_name.(cs)
    core_idxs = [findfirst(isequal(s), AtomicLevels.noble_gases) for s in cores]
    any(isnothing, core_idxs) && return nothing
    return AtomicLevels.noble_gases[minimum(core_idxs)]
end

principal(o::Orbital) = o.n
orbangmom(o::Orbital) = o.ℓ

function multireference(cs::Vector{Configuration{O}}; jvalues=nothing) where O <: Orbital
    maxℓ = maximum(orbangmom(o) for c in cs for (o, _, _) in c)
    maxorbs = [maximum(o for c in cs for (o, _, _) in c if orbangmom(o) == ℓ) for ℓ = 0:maxℓ]
    CSFDefinition(CSFDefConfiguration.(cs), maxorbs, 0, jvalues)
end

allorbs(n::Integer) = [Orbital(n, ℓ) for ℓ = 0:n-1]

function csfdefs_vv(cs, n)
    @assert n >= maximum(principal(o) for c in cs for (o, _, _) in c)
    confs = map(cs) do c
        minelecs = map(c) do (_, nelec, state)
            if (state == :closed) || (state == :inactive)
                nelec
            else
                @assert state == :open
                max(0, nelec - 2)
            end
        end
        CSFDefConfiguration(c, minelecs)
    end
    CSFDefinition(confs, allorbs(n), 2)
end

function csfdefs_cv(cs, min_n, max_n)
    @assert length(cs) >= 1
    @assert max_n >= maximum(principal(o) for c in cs for (o, _, _) in c)
    core_orbitals = unique(o for c in cs for (o, _, _) in core(c))
    @assert all(o ∈ core(c) for o in core_orbitals for c in cs)
    filter!(o -> principal(o) >= min_n, core_orbitals)
    cds = CSFDefinition[]
    for core_orbital in core_orbitals
        cdcs = CSFDefConfiguration[]
        for c in cs
            minelecs = map(c) do (o, nelec, state)
                if o == core_orbital
                    nelec - 1
                elseif (state == :closed) || (state == :inactive)
                    nelec
                else
                    @assert state == :open
                    max(0, nelec - 2)
                end
            end
            push!(cdcs, CSFDefConfiguration(c, minelecs))
        end
        push!(cds, CSFDefinition(cdcs, allorbs(max_n), 2))
    end
    return cds
end

function csfdefs_cc(cs, min_n, max_n)
    @assert length(cs) >= 1
    @assert max_n >= maximum(principal(o) for c in cs for (o, _, _) in c)
    core_orbitals = unique(o for c in cs for (o, _, _) in core(c))
    @assert all(o ∈ core(c) for o in core_orbitals for c in cs)
    filter!(o -> principal(o) >= min_n, core_orbitals)
    confs = map(cs) do c
        minelecs = map(c) do (o, nelec, state)
            if (state == :open) || (o in core_orbitals)
                max(0, nelec - 2)
            end
        end
        CSFDefConfiguration(c, minelecs)
    end
    CSFDefinition(confs, allorbs(max_n), 2)
end

function write_rcsfinput(io::IO, csfdefs; orbital_order="*", jvalues=nothing, core=nothing)
    # jvalues is allowed to be nothing, but otherwise we must be able to iterate over it
    if isnothing(jvalues)
        # If jvalues === nothing, then jvalues has to be provided for each CSF definition
        # separately
        if any(isnothing(csfdef.jvalues) for csfdef in csfdefs)
            throw(ArgumentError("some csfdefs are missing jvalues and jvalues kwarg not provided"))
        end
    else
        length(jvalues) == 2 || throw(ArgumentError("there must be two jvalues (passed: $jvalues)"))
        j1, j2 = try
            j1, j2 = jvalues
        catch e
            @error "Error occurred while unpacking jvalues" e
            throw(ArgumentError("unable to unpack jvalues (passed: $jvalues)"))
        end
        isinteger(2*j1) && isinteger(2*j2) || throw(ArgumentError("jvalues must be (half-)integers (passed: $jvalues)"))
        j1 >= 0 && j2 >= 0 || throw(ArgumentError("jvalues must be positive (passed: $jvalues)"))
        j2 >= j1 || throw(ArgumentError("first jvalues in the wrong order (passed: $jvalues)"))
    end

    core_idx = isnothing(core) ? 0 : findfirst(isequal(core), noble_cores)

    println(io, "$(orbital_order)")
    println(io, "$(core_idx)")
    for (i, csfdef) in enumerate(csfdefs)
        for c in csfdef.configurations
            for (o, maxelec, minelec) in c
                s = if maxelec == minelec
                    "i"
                elseif minelec == 0
                    "*"
                else
                    string(minelec)
                end
                print(io, "$(o)($maxelec,$s)")
            end
            println(io)
        end
        println(io) # newline indicating end of list of confs
        println(io, join(string.(csfdef.maxorbitals), ',')) # max orbitals?

        j1, j2 = if isnothing(jvalues)
            csfdef.jvalues
        else
            if !isnothing(csfdef.jvalues)
                @warn "jvalues in CSFDefinition #$i overriden by keyword argument" csfdef.jvalues jvalues
            end
            jvalues
        end
        min2J, max2J = convert(Int, 2*j1), convert(Int, 2*j2)
        println(io, "$(min2J), $(max2J)")

        println(io, csfdef.nsubstitutions) # number of excitations?
        println(io, i == length(csfdefs) ? "n" : "y") # another list?
    end
end

end
