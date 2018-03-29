module GRASP
using DocStringExtensions
import Humanize

export nelectrons, nexcitations

#
# Includes and imports related to compiled Fortran code.
# ------------------------------------------------------------------------------------------

let grasp_path_file = joinpath(dirname(@__FILE__), "../deps/grasp-path.jl")
    isfile(grasp_path_file) || error("deps/grasp-path.jl does not exist. Run `Pkg.build(\"GRASP\")`.")
    include(grasp_path_file)
end

const libgrasp_so = joinpath(dirname(@__FILE__), "..", "deps", "libgrasp.so")
isfile(libgrasp_so) || error("$(libgrasp_so) does not exist. Run `Pkg.build(\"GRASP\")`.")

#
# Additional includes and imports.
# ------------------------------------------------------------------------------------------

include("Symmetries.jl")
using .Symmetries

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

include("rcsfs.jl")
include("csfs.jl")
include("rmix.jl")

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
h2k_humanize(x) = Humanize.digitsep(round(Int, hartree2kayser(x)))

end # module
