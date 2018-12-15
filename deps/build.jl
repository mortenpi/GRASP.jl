if !("GRASP" in keys(ENV))
    error("\$GRASP environment variable not defined -- unable to determine the location of GRASP")
end

open("grasp-path.jl", "w") do io
    write(io, "const grasp = \"$(ENV["GRASP"])\"")
end
include("grasp-path.jl")

# Creating basic GRASP input/output files, using GRASP
# ./build-testfiles.sh cds into ../test/graspfiles
println("Running ./build-testfiles.sh")
run(`./build-testfiles.sh`)

# Build the Fortran shared library
run(`make`)
