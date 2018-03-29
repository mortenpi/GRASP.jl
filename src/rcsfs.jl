# Relativistic CSLs -- to parse rcsf.inp/out files.

struct FilledOrbital
    n :: Int
    kappa :: Int
    total2J :: Int
    nelectrons :: Int

    function FilledOrbital(n, kappa, total2J, nelectrons)
        @assert n >= 1
        @assert kappa != 0
        @assert kappa2l(kappa) < n
        new(n, kappa, total2J, nelectrons)
    end
end

struct CSF
    orbitals :: Vector{FilledOrbital}
    coupled2Js :: Vector{Int}
    total2J :: Int
    "State with positive or negative parity? `true` = positive."
    parity :: Bool
    CSF(total2J, parity, orbs, coupled2Js) = new(orbs, coupled2Js, total2J, parity)
end

struct CSFBlock
    csfs :: Vector{CSF}
    angularsym :: AngularSymmetry
    function CSFBlock(csfs::Vector{CSF})
        parity, total2J = first(csfs).parity, first(csfs).total2J
        @assert all(csf -> csf.parity == parity, csfs)
        @assert all(csf -> csf.total2J == total2J, csfs)
        new(csfs, AngularSymmetry(total2J, parity))
    end
end

total2J(cb::CSFBlock) = cb.angularsym.angmom.twoj
parity(cb::CSFBlock) = parity(cb.angularsym)

function Base.push!(csfb::CSFBlock, csf::CSF)
    @assert csf.total2J == total2J(csfb)
    push!(csfb.csfs, csf)
end

nelectrons(csf :: CSF) = sum(orb.nelectrons for orb in csf.orbitals)

"""
    nexcitations(csf_from::CSF, csf_to::CSF)

Calculates the excitations between relativistic orbitals. I.e. difference of occupation
number between `nl-` to `nl` counts as an excitation.
"""
function nexcitations(csf_from::CSF, csf_to::CSF)
    @assert nelectrons(csf_from) == nelectrons(csf_to)
    let from_orbs = [(orb.n, orb.kappa) for orb in csf_from.orbitals],
        to_orbs = [(orb.n, orb.kappa) for orb in csf_to.orbitals]
        @assert length(from_orbs) == length(unique(from_orbs))
        @assert length(to_orbs) == length(unique(to_orbs))
    end
    #from_nelecs = [orb.nelectrons for orb in csf_from.orbitals]

    to_unchecked = fill(true, length(csf_to.orbitals))

    nexcitations = 0
    nremoved = 0 # for sanity checking only

    # Let's count all the electrons that have been excited into some of
    # the occupied orbitals of csf_from.
    for orb in csf_from.orbitals
        toids = find(o -> o.n == orb.n && o.kappa == orb.kappa, csf_to.orbitals)
        to_nelectrons = if length(toids) == 0
            0
        elseif length(toids) == 1
            # Mark this csf_to orbital checked.
            toid = first(toids)
            @assert to_unchecked[toid]
            to_unchecked[toid] = false
            csf_to.orbitals[toid].nelectrons
        else
            error("Too many matches ($(length(toids))) among csf_to.orbitals.")
        end

        diff = to_nelectrons - orb.nelectrons
        # If diff < 0, electrons have been excited from this orbital
        (diff < 0) && (nremoved += abs(diff)) # for sanity checking
        nexcitations += (diff > 0) ? diff : 0
    end

    # Now let's loop over all the remaning csf_to orbitals
    for (i, orb) in enumerate(csf_to.orbitals)
        to_unchecked[i] || continue
        nexcitations += orb.nelectrons
    end

    @assert nexcitations == nremoved # sanity check
    return nexcitations
end

function parse_rcsf(filename)
    open(filename, "r") do io
        # First line should be "Core subshells:"
        line = readline(io); @assert strip(line) == "Core subshells:"
        line_cores = strip(readline(io), ['\n', '\r'])
        line = readline(io); @assert strip(line) == "Peel subshells:"
        line_peels = readline(io)
        line = readline(io); @assert strip(line) == "CSF(s):"

        core_orbs = parse_cores(line_cores)

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

            _total2J, _parity, orbs, coupled2Js = parse_csflines(line1, line2, line3)
            csf = CSF(
                _total2J,
                _parity,
                vcat(core_orbs, orbs),
                vcat(zeros(Int, length(core_orbs)), coupled2Js)
            )

            if isempty(csfblocks) || total2J(last(csfblocks)) != _total2J || parity(last(csfblocks)).even_parity != _parity
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

    orbs = FilledOrbital[]
    coupled2Js = Int[]
    for i = 1:div(length(line1), 9)
        orb = line1[9*(i-1)+1:9*i]
        @assert orb[6]=='(' && orb[9]==')'
        n = parse(Int, orb[1:3])
        kappa = parse_j(orb[4:5])
        nelec = parse(Int, orb[7:8])

        # pick out coupled angular momentum from the second line:
        orb2J = (length(line2) >= 9*i) ? parse2J(line2[9*(i-1)+1:9*i]) : 0

        # Pick the J-coupling from between the orbitals (line3).
        # The items in that line are shifted by three characters for some reason.
        if i > 1
            # TODO: Document this hack
            c2J_idx_first = 9*(i-1)+4
            c2J_idx_last = min(9*i+3, length(rstrip(line3))-1)
            # TODO: I used the following subrange at one point to fix one bug it appears it
            # created another one:
            #
            #     c2J_idx_first = 9*i+1
            #     c2J_idx_last  = min(9*(i+1), length(rstrip(line3))-1)
            #
            c2J_string = line3[c2J_idx_first:c2J_idx_last]
            c2J = try
                parse2J(c2J_string)
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
            push!(coupled2Js, c2J)
        else
            @assert strip(line3[9*(i-1)+2:9*i+1]) |> isempty
            # TODO: following related to the TODO above
            #@assert strip(line3[9*i+1:9*(i+1)]) |> isempty
        end

        push!(orbs, FilledOrbital(n, kappa, orb2J, nelec))
    end

    # Get the total angular momentum and parity
    parity = parse_parity(last(strip(line3)))
    total2J = last(coupled2Js)

    total2J, parity, orbs, coupled2Js
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


@deprecate parse2J(s) Int(parse(AngularMomentum, s).twoj)
@deprecate parse_parity(c::Char) Parity(c).even_parity
@deprecate parse_parity(s::AbstractString) parse(Parity, s).even_parity

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
    # Assuming NNNLL
    @show length(line)
    #@assert length(line) % 5 == 0

    orbstrings = split(line)
    orbs = FilledOrbital[]
    for os in orbstrings
        n, kappa = parse_orbital(os)
        push!(orbs, FilledOrbital(n, kappa, 0, 2*abs(kappa)))
    end
    orbs
end

# Base.print methods for the structs
function Base.print(io::IO, orb::FilledOrbital)
    write(io, string(orb.n), kappa2rso(orb.kappa))
    write(io, '(', string(orb.nelectrons), '~', string_2J(orb.total2J), ')')
    nothing
end

function Base.print(io::IO, csf::CSF)
    first = true
    for orb in csf.orbitals
        (! first) && write(io, ' ')
        print(io, orb)
        first = false
    end
    write(io, " ~ ", string_2J(csf.total2J))
    write(io, csf.parity ? '+' : '-')
    nothing
end

string_2J(twoJ::Integer) = (twoJ%2==0) ? "$(div(twoJ, 2))" : "$(twoJ)/2"
