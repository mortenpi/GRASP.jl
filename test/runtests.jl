using GRASP
using Base.Test

@testset "GRASP.jl" begin

include("csfs.jl")

@testset "libgrasp" begin
    @test isfile("test-libgrasp")
    @test isfile("graspfiles/rmix.out")

    # This checks that the library has been successfully linked against an executable.
    run(`./test-libgrasp graspfiles/rmix.out`)

    @test isa(GRASP.read_rmix("graspfiles/rmix.out"), GRASP.MixingFile)
end

end
