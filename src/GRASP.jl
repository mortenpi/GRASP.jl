module GRASP

let grasp_path_file = joinpath(dirname(@__FILE__), "../deps/grasp-path.jl")
    isfile(grasp_path_file) || error("deps/grasp-path.jl does not exist. Run `Pkg.build(\"GRASP\")`.")
    include(grasp_path_file)
end

include("csfs.jl")

end # module
