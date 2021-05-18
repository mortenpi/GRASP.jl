using Test
using GRASP

@testset "GRASP.jl" begin

@eval module SymmetriesTests include("symmetries.jl") end
@eval module CSFSTests include("Configurations.jl") end
@eval module RCSFSTests include("rcsfs.jl") end

@testset "libgrasp" begin
    # Make sure that all files exist
    @test isfile(joinpath(@__DIR__, "grasp/mixing/rmix.out"))

    @test_throws ErrorException GRASP.read_rmix("invalid-file.rmix")

    @test isa(
        GRASP.read_rmix(joinpath(@__DIR__, "grasp/mixing/rmix.out")),
        GRASP.MixingFile
    )

    @test isa(
        GRASP.read_rwfn(joinpath(@__DIR__, "grasp/mixing/rwfn.out")),
        Vector{GRASP.RWFNOrbital}
    )

    # make sure that existing, but invalid files throw gracefully
    @info "Expected error output:"
    @test_throws ErrorException GRASP.read_rmix(joinpath(@__DIR__, "grasp/mixing/rwfn.out"))
    @test_throws ErrorException GRASP.read_rwfn(joinpath(@__DIR__, "grasp/mixing/rmix.out"))
    @info "End of expected error output."

    include("binaries/binaries.jl")
end

@testset "other" begin
    @test GRASP.hartree2kayser(1.0) == 219474.63137
    @test GRASP.h2k_humanize(1.0) == "219,475"
end

end
