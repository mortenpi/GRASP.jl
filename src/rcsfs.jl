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
    total2J :: Int
    "State with positive or negative parity? `true` = positive."
    parity :: Bool
    CSF(total2J, parity, orbs) = new(orbs, total2J, parity)
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
        csfs = CSF[]
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

            #@show blockid, csfid
            #println(line1); println(line2); println(line3)
            total2J, parity, orbs = parse_csflines(line1, line2, line3)
            csf = CSF(total2J, parity, vcat(core_orbs, orbs))
            push!(csfs, csf)
        end
        return csfs
    end
end

function parse_csflines(line1, line2, line3)
    # Assuming that the CSF line consists of NNNLL(EE) blocks, each 9 characters long.
    @assert length(line1) % 9 == 0

    orbs = FilledOrbital[]
    for i = 1:div(length(line1), 9)
        orb = line1[9*(i-1)+1:9*i]
        @assert orb[6]=='(' && orb[9]==')'
        n = parse(Int, orb[1:3])
        kappa = parse_j(orb[4:5])
        nelec = parse(Int, orb[7:8])

        # pick out coupled angular momentum from the second line:
        orb2J = (length(line2) >= 9*i) ? parse2J(line2[9*(i-1)+1:9*i]) : 0

        push!(orbs, FilledOrbital(n, kappa, orb2J, nelec))
    end

    # Get the total angular momentum from line3

    total2J, parity = let x = strip(line3)
        parse2J(x[1:end-1]), parse_parity(x[end]) > 0
    end

    total2J, parity, orbs
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
    parse2J(s)

Parses the string `s` into the corresponding `2J`-value. String can either be a number or a
fraction of the form `<x>/2`/ An empty string (also one including just whitespace) parses
into 0.
"""
function parse2J(s)
    if in('/', s)
        J_numerator, J_denominator = split(s, '/'; limit=2)
        @assert parse(Int, J_denominator) == 2
        parse(Int, J_numerator)
    elseif !isempty(strip(s))
        2 * parse(Int, s)
    else
        0
    end
end

function parse_parity(c :: Char)
    if c == '+'
         1
    elseif c == '-'
        -1
    else
        error("Bad symbol $c")
    end
end

function parse_parity(s)
    c = strip(s)
    @assert length(c) == 1
    parse_parity(c[1])
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
