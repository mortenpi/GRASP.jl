"""
Types and methods related to parity.
"""
module Parities

export Parity, parity

"""
Represents the parity (even/odd) of a state.

Can be constructed from:

* `true` or `false` values, corresponding to even and odd, respectively.
* Characters `+` or `-`, corresponding to even and odd, respectively.
* The [`evenp`](@ref) and [`oddp`](@ref) constant instances.
"""
struct Parity
    even_parity :: Bool
end

"""
    even <: Parity

Represents even parity.
"""
const even = Parity(true)

"""
    odd <: Parity

Represents odd parity.
"""
const odd = Parity(false)

function Parity(c::Char)
    if c == '+'
        even
    elseif c == '-'
        odd
    else
        throw(ArgumentError("Invalid parity character '$c'. Must be '+' or '-'."))
    end
end

function Parity(s::AbstractString)
    if length(s) != 1
        throw(ArgumentError("Invalid argument `$s`. Must be a single character string."))
    end
    Parity(first(s))
end

Base.show(io::IO, p::Parity) = print(io, p.even_parity ? "+" : "-")

"""
    parity(object) -> Parity

Returns the parity (represented by an object of type [`Parity`](@ref)) of an object.
"""
function parity end

end
