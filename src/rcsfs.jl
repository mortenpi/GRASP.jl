# Relativistic CSLs -- to parse rcsf.inp/out files.

#
# type RelativisticOrbital
# ------------------------------------------------------------------------------------------
"""
    $(TYPEDEF)

Represents the label of a relativistic orbitals `nl(j/2)`. Internally it is
represented using the `n` and `kappa` values.

`isless` implements the conventional ordering on the orbitals.
"""
struct RelativisticOrbital
    n :: Int
    kappa :: Int

    function RelativisticOrbital(n, kappa)
        @assert n >= 1
        @assert kappa != 0
        @assert kappa2l(kappa) < n
        new(n, kappa)
    end
end

angularmomentum(ro::RelativisticOrbital) = AngularMomentum((2*abs(ro.kappa) - 1) // 2)

function Base.isless(ro1::RelativisticOrbital, ro2::RelativisticOrbital)
    if ro1.n != ro2.n
        return ro1.n < ro2.n
    elseif abs(ro1.kappa) != abs(ro2.kappa)
        return abs(ro1.kappa) < abs(ro2.kappa)
    else
        return ro1.kappa < ro2.kappa
    end
end

function Base.print(io::IO, ro::RelativisticOrbital)
    print(io, ro.n)
    write(io, kappa2rso(ro.kappa))
    nothing
end

maxelectrons(ro::RelativisticOrbital) = 2*abs(ro.kappa)

#
# type CSF
# ------------------------------------------------------------------------------------------
"""
    $(TYPEDEF)

Represents a configuration state function.

[`CSF`](@ref) is iterable, returning a tuple with the elements
`(orbital::RelativisticOrbital, n_electrons, orbital_coupling, csf_coupling)` on each
iteration.
"""
struct CSF
    angularsym :: AngularSymmetry
    orbitals :: Vector{RelativisticOrbital}
    occupations :: Vector{Int}
    orbcouplings :: Vector{AngularMomentum}
    csfcouplings :: Vector{AngularMomentum}

    function CSF(
            orbitals::Vector{RelativisticOrbital},
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

Base.start(csf::CSF) = 0
Base.next(csf::CSF, s) = (csf.orbitals[s+1], csf.occupations[s+1], csf.orbcouplings[s+1], csf.csfcouplings[s+1]), s + 1
Base.done(csf::CSF, s) = (s >= length(csf))
Base.length(csf::CSF) = length(csf.orbitals)

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
    let from_orbs = [(i, orb.n, orb.kappa) for (i, orb) in enumerate(csf_from.orbitals)],
        to_orbs = [(i, orb.n, orb.kappa) for (i, orb) in enumerate(csf_to.orbitals)]
        @assert length(from_orbs) == length(unique(from_orbs))
        @assert length(to_orbs) == length(unique(to_orbs))
    end

    to_unchecked = fill(true, length(csf_to.orbitals))

    nexcitations = 0
    nremoved = 0 # for sanity checking only

    # Let's count all the electrons that have been excited into some of
    # the occupied orbitals of csf_from.
    for (i, orb) in enumerate(csf_from.orbitals)
        toids = find(o -> o.n == orb.n && o.kappa == orb.kappa, csf_to.orbitals)
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
        write(io, '(', string(nelec), '|', string(orbcoupling), '|', string(orbcoupling), ')')
        first = false
    end
    write(io, " | ", string(csf.angularsym))
    nothing
end


#
# type CSFBlock
# ------------------------------------------------------------------------------------------

struct CSFBlock
    csfs :: Vector{CSF}
    angularsym :: AngularSymmetry

    function CSFBlock(csfs::Vector{CSF})
        angsym = first(csfs).angularsym
        @assert all(csf -> csf.angularsym == angsym, csfs)
        new(csfs, angsym)
    end
end

Base.length(cb::CSFBlock) = length(cb.csfs)

parity(cb::CSFBlock) = parity(cb.angularsym)
angularmomentum(cb::CSFBlock) = angularmomentum(cb.angularsym)

function Base.push!(csfb::CSFBlock, csf::CSF)
    @assert csf.angularsym == csfb.angularsym
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
        core_occupations = map(maxelectrons, core_orbitals)

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
            csf = CSF(
                vcat(core_orbitals, orbitals),
                vcat(core_occupations, noccupations),
                vcat(core_couplings, orbcouplings),
                vcat(core_couplings, csfcouplings),
                parity_,
            )

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

    orbs = RelativisticOrbital[]
    orbs_nelectrons = Int[]
    orbs_orbcouplings = AngularMomentum[]
    orbs_csfcouplings = AngularMomentum[]
    norbitals = div(length(line1), 9) # number of orbitals in this group
    for i = 1:norbitals
        orb = line1[9*(i-1)+1:9*i]
        @assert orb[6]=='(' && orb[9]==')'
        n = parse(Int, orb[1:3])
        kappa = parse_j(orb[4:5])
        nelec = parse(Int, orb[7:8])

        # pick out coupled angular momentum from the second line:
        angmom = (length(line2) >= 9*i) ? parse(AngularMomentum, line2[9*(i-1)+1:9*i]) : AngularMomentum(0)

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
            parse(AngularMomentum, c2J_string)
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

        push!(orbs, RelativisticOrbital(n, kappa))
        push!(orbs_nelectrons, nelec)
        push!(orbs_orbcouplings, angmom)
    end
    @assert length(orbs_csfcouplings) == norbitals

    # Get the total angular momentum and parity
    parity = Parity(last(strip(line3)))
    parity, orbs, orbs_nelectrons, orbs_orbcouplings, orbs_csfcouplings
end


# Helper functions for parsing
# TODO: naming here requires more work..

parse_l(lstring) = findfirst(ls -> ls == strip(lstring), SPECTROSCOPIC_NAMES) - 1

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

"""
Returns `(n, kappa)`.
"""
function parse_orbital(s)
    s = strip(s)
    if endswith(s,"-")
        l = parse_l(s[end-1:end-1])
        @assert l != 0
        parse(Int, s[1:end-2]), l
    else
        parse(Int, s[1:end-1]), -parse_l(s[end:end])-1
    end
end

function parse_cores(line)
    orbstrings = split(line)
    orbs = RelativisticOrbital[]
    for os in orbstrings
        n, kappa = parse_orbital(os)
        push!(orbs, RelativisticOrbital(n, kappa))
    end
    orbs
end
