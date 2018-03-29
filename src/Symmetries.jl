"""
Types and methods related to parity.
"""
module Symmetries

export Parity, parity
export AngularMomentum, angularmomentum
export AngularSymmetry

#
# Parity
# ------------------------------------------------------------------------------------------

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

#
# Angular momentum
# ------------------------------------------------------------------------------------------

"""
Represent integer and half-integer angular momentum values.

Defined through the `2J` value.
"""
struct AngularMomentum
    twoj :: UInt
    function AngularMomentum(twoj :: Integer)
        if twoj < 0
            throw(ArgumentError("Angular momentum needs to be positive."))
        end
        new(twoj)
    end
end

Base.show(io::IO, am::AngularMomentum) = print(io, iseven(am.twoj) ? string(div(am.twoj, 2)) : "$(am.twoj)/2")

"""
    angularmomentum(object) -> AngularMomentum

Returns the angular momentum of an object, represented by an [`AngularMomentum`](@ref)
object.
"""
function angularmomentum end

#
# Angular symmetry (angular momentum + parity)
# ------------------------------------------------------------------------------------------

"""
Represents the angular symmetry (parity and angular momentum) of a state.
"""
struct AngularSymmetry
    angmom :: AngularMomentum
    parity :: Parity
    AngularSymmetry(angmom::AngularMomentum, parity::Parity) = new(angmom, parity)
end
AngularSymmetry(angmom, parity::Parity) = AngularSymmetry(AngularMomentum(angmom), parity)
AngularSymmetry(angmom, parity) = AngularSymmetry(angmom, Parity(parity))

angularmomentum(as::AngularSymmetry) = as.angmom
parity(as::AngularSymmetry) = as.parity

function Base.show(io::IO, as::AngularSymmetry)
    show(io, as.angmom)
    show(io, as.parity)
end

end
