# Relativistic CSFs

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
    orbitals :: Vector{RelativisticOrbital{Int}}
    occupations :: Vector{Int}
    orbcouplings :: Vector{AngularMomentum}
    csfcouplings :: Vector{AngularMomentum}

    function CSF(
            orbitals::Vector{<:RelativisticOrbital},
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

function kappa2rso(kappa :: Integer)
    lstr = specname(kappa_to_l(kappa))
    (kappa < 0) ? lstr : lstr*"-"
end
