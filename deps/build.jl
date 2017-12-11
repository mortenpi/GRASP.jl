if !("GRASP2K" in keys(ENV))
    error("\$GRASP2K environment variable not defined -- unable to determine the location of GRASP2K")
end

open("grasp-path.jl", "w") do io
    write(io, "const grasp2k = \"$(ENV["GRASP2K"])\"")
end
include("grasp-path.jl")

# Creating basic GRASP input/output files, using GRASP
# ./build-testfiles.sh cds into ../test/graspfiles
println("Running ./build-testfiles.sh")
run(`./build-testfiles.sh`)

# Build the Fortran shared library
run(`make`)
