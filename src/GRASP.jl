module GRASP
using DocStringExtensions

export angularmomentum, parity # from module Symmetries
export nelectrons, nexcitations, maxelectrons

import AtomicLevels
using AtomicLevels: Orbital, RelativisticOrbital, degeneracy
using HalfIntegers: HalfInt

#
# Includes and imports related to compiled Fortran code.
# ------------------------------------------------------------------------------------------

let grasp_path_file = joinpath(dirname(@__FILE__), "../deps/grasp-path.jl")
    isfile(grasp_path_file) || error("deps/grasp-path.jl does not exist. Run `Pkg.build(\"GRASP\")`.")
    include(grasp_path_file)
end

const libgraspjl_so = joinpath(dirname(@__FILE__), "..", "deps", "libgraspjl.so")
isfile(libgraspjl_so) || error("$(libgrasp_so) does not exist. Run `Pkg.build(\"GRASP\")`.")

const libgrasp_so = joinpath(dirname(@__FILE__), "..", "deps", "libgrasp.so")
isfile(libgrasp_so) || error("$(libgrasp_so) does not exist. Run `Pkg.build(\"GRASP\")`.")

#
# Additional includes and imports.
# ------------------------------------------------------------------------------------------

include("Symmetries.jl")
import .Symmetries: Parity, AngularMomentum, AngularSymmetry
import .Symmetries: parity, angularmomentum, absdiff

const SPECTROSCOPIC_NAMES = "s p d f g h i k l m n o q r t u v" |> split
function specname(l::Integer)
    @assert l >= 0
    SPECTROSCOPIC_NAMES[l+1]
end

"""
    nelectrons(obj)

Calculate the number of electrons in `obj`.
"""
function nelectrons end

"""
    nexcitations(from, to)

Calculate the number of excitations needed to go from `from` to `to`.
"""
function nexcitations end

"""
    maxelectrons(obj)

Return the maximum number of electrons.
"""
function maxelectrons end

include("rcsfs.jl")
include("Configurations.jl")
include("rmix.jl")
include("rwfn.jl")
include("rcsfparser.jl")

# Fixed Humanize.digitsep --- fixed for negative numbers
# from: https://github.com/IainNZ/Humanize.jl
function digitsep(value::Integer, sep = ",", k = 3)
    isnegative = value < 0
    value = string(abs(value))
    n = length(value)
    starts = reverse(collect(n:-k:1))
    groups = [value[max(x-k+1, 1):x] for x in starts]
    return (isnegative ? "-" : "") * join(groups, sep)
end

"""
    $(SIGNATURES)

Converts the energy value `x` from Hartrees to Kaysers (aka 1/cm),
using the conversion ``1~\\textrm{Ha} = 219474.63137~\\textrm{cm}^{-1}``.

# Examples

```jldoctest
julia> GRASP.hartree2kayser(1.0)
219474.63137
```
"""
hartree2kayser(x) = 219474.63137*x

"""
    $(SIGNATURES)

Like [`hartree2kayser`](@ref), but produces a humanized string with thousands separators.

# Examples

```jldoctest
julia> GRASP.h2k_humanize(1.0)
"219,475"
```
"""
h2k_humanize(x) = digitsep(round(Int, hartree2kayser(x)))

end # module
