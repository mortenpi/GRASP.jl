if !("GRASP2K" in keys(ENV))
    error("\$GRASP2K environment variable not defined -- unable to determine the location of GRASP2K")
end

open("grasp-path.jl", "w") do io
    write(io, "const grasp2k = \"$(ENV["GRASP2K"])\"")
end
