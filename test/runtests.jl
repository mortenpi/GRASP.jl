using GRASP
using Base.Test

@testset "GRASP.jl" begin

include("csfs.jl")
include("rcsfs.jl")

@testset "libgrasp" begin
    @test isfile("test-libgrasp")
    @test isfile("graspfiles/rmix.out")

    # This checks that the library has been successfully linked against an executable.
    run(`./test-libgrasp graspfiles/rmix.out`)

    @test isa(GRASP.read_rmix("graspfiles/rmix.out"), GRASP.MixingFile)
end

@testset "other" begin
    @test GRASP.hartree2kayser(1.0) == 219474.63137
    @test GRASP.h2k_humanize(1.0) == "219,475"
end

end
