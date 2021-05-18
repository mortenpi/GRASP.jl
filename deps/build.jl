"GRASP" in keys(ENV) || error("\$GRASP environment variable not defined -- unable to determine the location of GRASP")
isabspath(ENV["GRASP"]) || error("\$GRASP not set to an absolute path (GRASP=$(ENV["GRASP"]))")
isdir(ENV["GRASP"]) || error("\$GRASP not pointing to a directory (GRASP=$(ENV["GRASP"]))")

open("grasp-path.jl", "w") do io
    grasp_path = abspath(ENV["GRASP"])
    write(io, "const grasp = \"$(grasp_path)\"")
end
include("grasp-path.jl")

# Creating basic GRASP input/output files, using GRASP
# ./build-testfiles.sh cds into ../test/graspfiles
println("Running ./build-testfiles.sh")
run(`./build-testfiles.sh`)

# Build the Fortran shared library
run(`make`)
