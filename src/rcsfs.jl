# Relativistic CSLs -- to parse rcsf.inp/out files.

"""
    kappa_to_l(κ)

Calculate the l quantum number corresponding to the κ quantum number.

Note: κ and l values are always integers.
"""
function kappa_to_l(kappa::Integer)
    kappa == zero(kappa) && throw(ArgumentError("κ can not be zero"))
    (kappa < 0) ? -(kappa+1) : kappa
end

function kappa(orb::Orbital{I,Rational{I}}) where {I}
    (orb.j < orb.ℓ ? 1 : -1) * I(orb.j + 1//2)
end

#
# type CSF
# ------------------------------------------------------------------------------------------
"""
    $(TYPEDEF)

Represents a configuration state function.

[`CSF`](@ref) is iterable, returning a tuple with the elements
`(orbital::Orbital, n_electrons, orbital_coupling, csf_coupling)` on each
iteration.
"""
struct CSF
    angularsym :: AngularSymmetry
    orbitals :: Vector{Orbital{Int,Rational{Int}}}
    occupations :: Vector{Int}
    orbcouplings :: Vector{AngularMomentum}
    csfcouplings :: Vector{AngularMomentum}

    function CSF(
            orbitals::Vector{Orbital{Int,Rational{Int}}},
            electrons::Vector{Int},
            couplings_orbitals::Vector{AngularMomentum},
            couplings_csf::Vector{AngularMomentum},
            parity,
        )
        let norbitals = length(orbitals)
            @assert length(electrons) == norbitals
            @assert length(couplings_orbitals) == norbitals
            @assert length(couplings_csf) == norbitals
        end
        @assert issorted(orbitals)
        @assert all(x -> x!=0, electrons)
        new(
            AngularSymmetry(last(couplings_csf), parity),
            orbitals, electrons, couplings_orbitals, couplings_csf,
        )
    end
end

import Base: ==
function ==(csf1::CSF, csf2::CSF)
    (csf1.angularsym == csf2.angularsym)     &&
    (csf1.orbitals == csf2.orbitals)         &&
    (csf1.occupations == csf2.occupations)   &&
    (csf1.orbcouplings == csf2.orbcouplings) &&
    (csf1.csfcouplings == csf2.csfcouplings)
end

function Base.iterate(csf::CSF, s=1)
    s > length(csf) && return nothing
    item = (csf.orbitals[s], csf.occupations[s], csf.orbcouplings[s], csf.csfcouplings[s])
    return item, s + 1
end

Base.length(csf::CSF) = length(csf.orbitals)

angularmomentum(orb::Orbital) = AngularMomentum(orb.j)

parity(csf::CSF) = parity(csf.angularsym)
angularmomentum(csf::CSF) = angularmomentum(csf.angularsym)
nelectrons(csf :: CSF) = sum(csf.occupations)

"""
    nexcitations(csf_from::CSF, csf_to::CSF)

Calculates the excitations between relativistic orbitals. I.e. difference of occupation
number between `nl-` to `nl` counts as an excitation.
"""
function nexcitations(csf_from::CSF, csf_to::CSF)
    @assert nelectrons(csf_from) == nelectrons(csf_to)
    let from_orbs = [(i, orb.n, kappa(orb)) for (i, orb) in enumerate(csf_from.orbitals)],
        to_orbs = [(i, orb.n, kappa(orb)) for (i, orb) in enumerate(csf_to.orbitals)]
        @assert length(from_orbs) == length(unique(from_orbs))
        @assert length(to_orbs) == length(unique(to_orbs))
    end

    to_unchecked = fill(true, length(csf_to.orbitals))

    nexcitations = 0
    nremoved = 0 # for sanity checking only

    # Let's count all the electrons that have been excited into some of
    # the occupied orbitals of csf_from.
    for (i, orb) in enumerate(csf_from.orbitals)
        toids = findall(o -> o.n == orb.n && kappa(o) == kappa(orb), csf_to.orbitals)
        to_nelectrons = if length(toids) == 0
            0
        elseif length(toids) == 1
            # Mark this csf_to orbital checked.
            toid = first(toids)
            @assert to_unchecked[toid]
            to_unchecked[toid] = false
            csf_to.occupations[toid]
        else
            error("Too many matches ($(length(toids))) among csf_to.orbitals.")
        end

        diff = to_nelectrons - csf_from.occupations[i]
        # If diff < 0, electrons have been excited from this orbital
        (diff < 0) && (nremoved += abs(diff)) # for sanity checking
        nexcitations += (diff > 0) ? diff : 0
    end

    # Now let's loop over all the remaning csf_to orbitals
    for (i, orb) in enumerate(csf_to.orbitals)
        to_unchecked[i] || continue
        nexcitations += csf_to.occupations[i]
    end

    @assert nexcitations == nremoved # sanity check
    return nexcitations
end

function Base.show(io::IO, csf::CSF)
    write(io, '"')
    print(io, csf)
    write(io, '"')
end

function Base.print(io::IO, csf::CSF)
    first = true
    for (orb, nelec, orbcoupling, csfcoupling) in csf
        !first && write(io, ' ')
        print(io, orb)
        write(io, '(', string(nelec), '|', string(orbcoupling), '|', string(csfcoupling), ')')
        first = false
    end
    write(io, " | ", string(csf.angularsym))
    nothing
end


#
# type CSFBlock
# ------------------------------------------------------------------------------------------

const AtLevCSF = AtomicLevels.CSF{Int,Rational{Int},Rational{Int}}
struct CSFBlock
    csfs :: Vector{AtLevCSF}
    angularsym :: AngularSymmetry

    function CSFBlock(csfs::Vector{AtLevCSF})
        csf = first(csfs)
        angsym = angularsym(csf)
        @assert all(csf -> angularsym(csf) == angsym, csfs)
        new(csfs, angsym)
    end
end

function angularsym(csf::AtomicLevels.CSF)
     AngularSymmetry(angularmomentum(csf), Parity(AtomicLevels.parity(csf.config)))
end

function Symmetries.Parity(x::Integer)
    if x == 1
        return Parity(true)
    elseif x == -1
        return Parity(false)
    else
        throw(ArgumentError("Invalid numeric parity value $x."))
    end
end

angularmomentum(csf::AtomicLevels.CSF) = AngularMomentum(last(csf.terms))

Base.length(cb::CSFBlock) = length(cb.csfs)

parity(cb::CSFBlock) = parity(cb.angularsym)
angularmomentum(cb::CSFBlock) = angularmomentum(cb.angularsym)

function Base.push!(csfb::CSFBlock, csf::AtLevCSF)
    @assert angularsym(csf) == csfb.angularsym
    push!(csfb.csfs, csf)
end

function parse_rcsf(filename)
    open(filename, "r") do io
        # First line should be "Core subshells:"
        line = readline(io); @assert strip(line) == "Core subshells:"
        line_cores = strip(readline(io), ['\n', '\r'])
        line = readline(io); @assert strip(line) == "Peel subshells:"
        line_peels = readline(io)
        line = readline(io); @assert strip(line) == "CSF(s):"

        core_orbitals = parse_cores(line_cores)
        core_couplings = AngularMomentum[0 for co in core_orbitals]
        core_occupations = map(degeneracy, core_orbitals)

        blockid, csfid = 1, 1
        csfblocks = CSFBlock[]
        while ! eof(io)
            line1 = readline(io)
            if startswith(line1, " *")
                blockid += 1
                csfid = 1
                line1 = readline(io)
            end
            line1 = strip(line1, ['\n', '\r'])
            line2 = strip(readline(io), ['\n', '\r'])
            line3 = strip(readline(io), ['\n', '\r'])

            parity_, orbitals, noccupations, orbcouplings, csfcouplings = parse_csflines(line1, line2, line3)

            @assert !any(isequal(0), noccupations) # we should never get 0-electron orbitals

            # Fix orbcouplings that are not explicitly written in the CSL file (represented
            # with nothings in the array). It is assumed that this is the case only if the
            # orbital is fully occupied.
            #
            # NOTE: we could, though, also omit values if there is only 1 electron or
            # maxelectrons(orb) - 1 electrons, this which case the angular momentum can
            # only have only one value too.
            for i = 1:length(orbitals)
                orbital, nelec = orbitals[i], noccupations[i]
                if orbcouplings[i] === nothing
                    nelec == degeneracy(orbital) || error("""
                    Unable to fix missing orbital coupling.
                      Orbital $i, orbital=$(repr(orbital)), nelec=$nelec
                    1: $(line1)
                    2: $(line2)
                    3: $(line3)
                    """)
                    orbcouplings[i] = 0
                elseif (nelec == 1 || nelec == degeneracy(orbital) - 1) && orbcouplings[i] != angularmomentum(orbital)
                    @warn "Bad orbital coupling" orbcouplings[i] nelec angularmomentum(orbital)
                elseif orbcouplings[i] > nelec * angularmomentum(orbital)
                    # If the orbital coupling is larger than (# particles) * (l of orbital),
                    # then that has to be an error.
                    @warn "Bad orbital coupling" orbcouplings[i] nelec angularmomentum(orbital)
                end
            end

            # Fix csfcouplings which are not explicitly written to the CSL file. This appears
            # to be the case for "trivial" zero couplings, if the previous CSF layer and
            # current orbital are both zeros.
            for i = 1:length(orbitals)
                oj = orbcouplings[i]
                cj = (i > 1) ? csfcouplings[i-1] : AngularMomentum(0)
                Δupper, Δlower = oj+cj, absdiff(oj,cj)
                if csfcouplings[i] === nothing
                    Δupper == Δlower || error("""
                    Unable to fix missing CSF coupling.
                      Orbital $i, orbital=$(repr(orbitals[i])), oj=$(repr(oj)), cj=$(repr(cj))
                    1: $(line1)
                    2: $(line2)
                    3: $(line3)
                    """)
                    csfcouplings[i] = Δlower
                elseif !(Δlower <= csfcouplings[i] <= Δupper)
                    @warn "Invalid csfcoupling value?" csfcouplings[i] Δupper Δlower
                end
            end

            config = AtomicLevels.Configuration(vcat(core_orbitals, orbitals), vcat(core_occupations, noccupations))
            subshell_terms = map(x -> convert(Rational{Int}, x),
                vcat(core_couplings, Vector{AngularMomentum}(orbcouplings)))
            terms = map(x -> convert(Rational{Int}, x),
                vcat(core_couplings, Vector{AngularMomentum}(csfcouplings)))
            csf = AtomicLevels.CSF(config, subshell_terms, terms)

            if isempty(csfblocks) || angularmomentum(last(csfblocks)) != last(csfcouplings) || parity(last(csfblocks)) != parity_
                push!(csfblocks, CSFBlock([csf]))
            else
                push!(last(csfblocks), csf)
            end
        end
        return csfblocks
    end
end

function parse_csflines(line1, line2, line3)
    # Assuming that the CSF line consists of NNNLL(EE) blocks, each 9 characters long.
    @assert length(line1) % 9 == 0

    orbs = Orbital{Int,Rational{Int}}[]
    orbs_nelectrons = Int[]
    orbs_orbcouplings = Union{AngularMomentum,Nothing}[]
    orbs_csfcouplings = Union{AngularMomentum,Nothing}[]
    norbitals = div(length(line1), 9) # number of orbitals in this group
    for i = 1:norbitals
        orb = line1[9*(i-1)+1:9*i]
        @assert orb[6]=='(' && orb[9]==')'
        orbital = AtomicLevels.orbital_from_string((orb[1:5]))
        # n = parse(Int, orb[1:3])
        # kappa = parse_j(orb[4:5])
        nelec = parse(Int, orb[7:8])

        # pick out coupled angular momentum from the second line:
        angmom =  if length(line2) < 9*i
            nothing
        else
            parse_angmom_string(line2[9*(i-1)+1:9*i])
        end

        # Pick the J-coupling from between the orbitals (line3).
        # The items in that line are shifted by three characters to the right for
        # some reason.. except for the last one, which defines the J^P values of
        # the CSF. That one is shifted by 1 character.. and then one after that
        # is the parity +/- symbol.
        c2J_idx_first, c2J_idx_last = if i < norbitals
            # So, for the non-last ones we assume a 9 character block that has been
            # shifted by 3 characters to the right.
            # +1 to go from 0-based to 1-based
            9*(i-1)+3+1, 9*i+3
        else # i == norbitals -- the last non-regular coupling
            9*(i-1)+3+1, 9*i+1
        end
        c2J_string = line3[c2J_idx_first:c2J_idx_last]
        coupled_angmom = try
            parse_angmom_string(c2J_string)
        catch e
            error("""
            Error parsing 2J value on line 3 (i=$i)
              range $(c2J_idx_first):$(c2J_idx_last) -> '$(c2J_string)'
            1: $(line1)
            2: $(line2)
            3: $(line3)
            $(' '^(c2J_idx_first+2))^$('-'^(c2J_idx_last-c2J_idx_first-1))^
            """)
        end
        push!(orbs_csfcouplings, coupled_angmom)

        push!(orbs, orbital)
        push!(orbs_nelectrons, nelec)
        push!(orbs_orbcouplings, angmom)
    end
    @assert length(orbs_csfcouplings) == norbitals

    # Get the total angular momentum and parity
    parity = Parity(last(strip(line3)))
    parity, orbs, orbs_nelectrons, orbs_orbcouplings, orbs_csfcouplings
end

function parse_angmom_string(s)
    length(strip(s)) == 0 && return nothing
    parse(AngularMomentum, s)
end

# Helper functions for parsing
# TODO: naming here requires more work..

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

kappa2l(kappa :: Integer) = (kappa < 0) ? -(kappa+1) : kappa

function kappa2rso(kappa :: Integer)
    lstr = specname(kappa2l(kappa))
    (kappa < 0) ? lstr : lstr*"-"
end

function parse_cores(line)
    orbstrings = split(line)
    orbs = Orbital{Int,Rational{Int}}[]
    for os in orbstrings
        push!(orbs, AtomicLevels.orbital_from_string(os))
    end
    orbs
end
