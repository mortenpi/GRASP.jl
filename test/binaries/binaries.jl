# This checks that the library has been successfully linked against an executable.
test_libgrasp_bin = joinpath(@__DIR__, "test-libgrasp")
test_libgrasp_arg = joinpath(@__DIR__, "../grasp/mixing/rmix.out")
@test success(`$(test_libgrasp_bin) $(test_libgrasp_arg)`)
