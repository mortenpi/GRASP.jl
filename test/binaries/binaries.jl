# This checks that the library has been successfully linked against an executable.
@testset "test-libgrasp" begin
    test_libgrasp_bin = joinpath(@__DIR__, "test-libgrasp")
    @test isfile(test_libgrasp_bin)

    test_libgrasp_arg = joinpath(@__DIR__, "../grasp/mixing/rmix.out")
    @test isfile(test_libgrasp_arg)
    @test success(`$(test_libgrasp_bin) $(test_libgrasp_arg)`)
end

@testset "test-readrwfn" begin
    test_readrwfn_bin = joinpath(@__DIR__, "test-readrwfn")
    @test isfile(test_readrwfn_bin)

    test_readrwfn_arg = joinpath(@__DIR__, "../grasp/mixing/rwfn.inp")
    @test isfile(test_readrwfn_arg)
    @test success(`$(test_readrwfn_bin) $(test_readrwfn_arg)`)

    test_readrwfn_arg = joinpath(@__DIR__, "../grasp/mixing/rwfn.out")
    @test isfile(test_readrwfn_arg)
    @test success(`$(test_readrwfn_bin) $(test_readrwfn_arg)`)
end
