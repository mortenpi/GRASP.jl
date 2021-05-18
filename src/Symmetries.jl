"""
Types and methods related to parity.
"""
module Symmetries

export Parity, parity
export AngularMomentum, angularmomentum
export AngularSymmetry
export absdiff

using AtomicLevels
using HalfIntegers: HalfInteger, HalfInt, half

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

function Parity(x::Integer)
    if x == 1
        return Parity(true)
    elseif x == -1
        return Parity(false)
    else
        throw(ArgumentError("Invalid numeric parity value $x."))
    end
end

Parity(x::AtomicLevels.Parity) = Parity(x.p)

function Base.parse(::Type{Parity}, s::AbstractString)
    s = strip(s)
    length(s) == 1 || throw(Meta.ParseError("Can't parse '$s' into Parity (bad length)."))
    try
        Parity(first(strip(s)))
    catch e
        isa(e, ArgumentError) && throw(Meta.ParseError(e.msg))
        rethrow()
    end
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

Use `Rational` values to define half-integer angular momenta.
"""
struct AngularMomentum
    twoj :: UInt
    function AngularMomentum(j::Union{Integer, Rational, HalfInt})
        if j < zero(j)
            throw(ArgumentError("Angular momentum can't be negative: $j"))
        end

        twoj = if denominator(j) == 1
            2*numerator(j)
        elseif denominator(j) == 2
            numerator(j)
        else
            throw(ArgumentError("Angular momentum needs to be integer or half-integer: $j"))
        end

        return new(twoj)
    end
end

Base.convert(::Type{AngularMomentum}, j::Union{Integer, Rational}) = AngularMomentum(j)
Base.convert(::Type{Rational}, am::AngularMomentum) = convert(Rational{Int}, am)
Base.convert(::Type{Rational{Int}}, am::AngularMomentum) = Rational(Int(am.twoj), 2)
Base.convert(::Type{HalfInteger}, am::AngularMomentum) = half(convert(Int, am.twoj))

import Base: +, *
+(a::AngularMomentum, b::AngularMomentum) = AngularMomentum(convert(Rational, a) + convert(Rational, b))
function *(a::T, b::AngularMomentum) where {T <: Integer}
    a < zero(T) && throw(ArgumentError("Can't multiply AngularMomentum with a negative number."))
    AngularMomentum(a * b.twoj // 2)
end
absdiff(a::AngularMomentum, b::AngularMomentum) = AngularMomentum(abs(convert(Rational, a) - convert(Rational, b)))

Base.isless(a::AngularMomentum, b::AngularMomentum) = (a.twoj < b.twoj)

"""
    parse(::Type{AngularMomentum}, ::AbstractString)

Parses the string `s` into the corresponding `2J`-value. String can either be a number or a
fraction of the form `<x>/2`/ An empty string (also one including just whitespace) parses
into 0.
"""
function Base.parse(::Type{AngularMomentum}, s::AbstractString)
    if in('/', s)
        J_numerator, J_denominator = split(s, '/'; limit=2)
        @assert parse(Int, J_denominator) == 2
        AngularMomentum(parse(Int, J_numerator) // 2)
    elseif !isempty(strip(s))
        AngularMomentum(parse(Int, s))
    else
        AngularMomentum(0)
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
    AngularSymmetry(angmom, parity::Parity) = new(angmom, parity)
end
AngularSymmetry(angmom, parity) = AngularSymmetry(angmom, Parity(parity))

angularmomentum(as::AngularSymmetry) = as.angmom
parity(as::AngularSymmetry) = as.parity

function Base.show(io::IO, as::AngularSymmetry)
    show(io, as.angmom)
    show(io, as.parity)
end

end
