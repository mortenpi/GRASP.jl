# This checks that the library has been successfully linked against an executable.
test_libgrasp_bin = joinpath(@__DIR__, "test-libgrasp")
test_libgrasp_arg = joinpath(@__DIR__, "../grasp/mixing/rmix.out")
@test success(`$(test_libgrasp_bin) $(test_libgrasp_arg)`)

test_readrwfn_bin = joinpath(@__DIR__, "test-readrwfn")
test_readrwfn_arg = joinpath(@__DIR__, "../grasp/mixing/rwfn.out")
@test success(`$(test_readrwfn_bin) $(test_readrwfn_arg)`)
