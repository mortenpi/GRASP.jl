module GRASP

let grasp_path_file = joinpath(dirname(@__FILE__), "../deps/grasp-path.jl")
    isfile(grasp_path_file) || error("deps/grasp-path.jl does not exist. Run `Pkg.build(\"GRASP\")`.")
    include(grasp_path_file)
end

const libgrasp_so = joinpath(dirname(@__FILE__), "..", "deps", "libgrasp.so")
isfile(libgrasp_so) || error("$(libgrasp_so) does not exist. Run `Pkg.build(\"GRASP\")`.")


const SPECTROSCOPIC_NAMES = "s p d f g h i k l m n o q r t u v" |> split
function specname(l::Integer)
    @assert l >= 0
    SPECTROSCOPIC_NAMES[l+1]
end


include("csfs.jl")
include("rcsfs.jl")
include("rmix.jl")

end # module
