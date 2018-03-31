using GRASP
using Base.Test

@testset "GRASP.jl" begin

@testset "Symmetries" begin include("symmetries.jl") end

include("csfs.jl")
include("rcsfs.jl")

@testset "libgrasp" begin
    @test isfile("test-libgrasp")
    @test isfile("grasp/mixing/rmix.out")

    # This checks that the library has been successfully linked against an executable.
    run(`./test-libgrasp grasp/mixing/rmix.out`)

    @test isa(GRASP.read_rmix("grasp/mixing/rmix.out"), GRASP.MixingFile)
end

@testset "other" begin
    @test GRASP.hartree2kayser(1.0) == 219474.63137
    @test GRASP.h2k_humanize(1.0) == "219,475"
end

end
