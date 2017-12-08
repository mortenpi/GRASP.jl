using Base.Test
import GRASP

import GRASP: CSFOrbital, CSFDefinition, CSFDefinitionList

@testset "csfs.jl" begin

@testset "specname" begin
    @test GRASP.specname(0) == "s"
    @test GRASP.specname(1) == "p"
    @test GRASP.specname(2) == "d"
    @test GRASP.specname(3) == "f"
end

@testset "CSFOrbital" begin
    @test CSFOrbital(3,0) == CSFOrbital(3,0)
    @test CSFOrbital(3,0) != CSFOrbital(3,1)
    @test CSFOrbital(3,0) != CSFOrbital(2,0)
    @test CSFOrbital(3,0) != CSFOrbital(2,1)
end

@testset "CSFDefinition" begin
    cdef1 = let cdef = CSFDefinition()
        for n=1:3, l=0:(n-1)
             push!(cdef, GRASP.CSFOrbital(n, l), 2*(2l+1), 0)
        end
        cdef
    end

    cdef2 = let cdef = CSFDefinition()
        for n=4:5, l=0:min(3, n-1)
             push!(cdef, GRASP.CSFOrbital(n, l), 2*(2l+1), 0)
        end
        cdef
    end

    @test string(cdef1) == "1s(2,i)2s(2,i)2p(6,i)3s(2,i)3p(6,i)3d(10,i)"
    @test string(cdef2) == "4s(2,i)4p(6,i)4d(10,i)4f(14,i)5s(2,i)5p(6,i)5d(10,i)5f(14,i)"

    @test string(cdef1 * cdef2) == "1s(2,i)2s(2,i)2p(6,i)3s(2,i)3p(6,i)3d(10,i)4s(2,i)4p(6,i)4d(10,i)4f(14,i)5s(2,i)5p(6,i)5d(10,i)5f(14,i)"
    # TODO: Test error for *?
end

end # @testset "csfs.jl"
