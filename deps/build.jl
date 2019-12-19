if !("GRASP" in keys(ENV))
    error("\$GRASP environment variable not defined -- unable to determine the location of GRASP")
end

_libgrasp_path = abspath(joinpath(ENV["GRASP"], "lib", "libgrasp-rci.so"))
if !isfile(_libgrasp_path)
    error("Unable to find libgrasp-rci.so (at $(_libgrasp_path))")
end

open("grasp-path.jl", "w") do io
    write(io, "const grasp = \"$(ENV["GRASP"])\"\n")
end
include("grasp-path.jl")

# Creating basic GRASP input/output files, using GRASP
# ./build-testfiles.sh cds into ../test/graspfiles
println("Running ./build-testfiles.sh")
run(`./build-testfiles.sh`)

# Build the Fortran shared library
run(`make`)
