using GRASP
using Base.Test

@testset "GRASP.jl" begin

@eval module SymmetriesTests include("symmetries.jl") end
@eval module CSFSTests include("csfs.jl") end
@eval module RCSFSTests include("rcsfs.jl") end

@testset "libgrasp" begin
    @test isfile(joinpath(@__DIR__, "test-libgrasp"))
    @test isfile(joinpath(@__DIR__, "grasp/mixing/rmix.out"))

    # This checks that the library has been successfully linked against an executable.
    test_libgrasp_script = joinpath(@__DIR__, "test-libgrasp")
    test_libgrasp_arg = joinpath(@__DIR__, "grasp/mixing/rmix.out")
    run(`$(test_libgrasp_script) $(test_libgrasp_arg)`)

    @test isa(
        GRASP.read_rmix(joinpath(@__DIR__, "grasp/mixing/rmix.out")),
        GRASP.MixingFile
    )
end

@testset "other" begin
    @test GRASP.hartree2kayser(1.0) == 219474.63137
    @test GRASP.h2k_humanize(1.0) == "219,475"
end

end
