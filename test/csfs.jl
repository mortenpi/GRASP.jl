using Base.Test
using GRASP
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

    # ordering
    orbitals = [CSFOrbital(n,l) for n=1:4 for l=0:(n-1)]
    for orb in orbitals
        @test orb == orb
    end

    for i=1:length(orbitals), j=(i+1):length(orbitals)
        a, b = orbitals[i], orbitals[j]
        @test a != b
        @test a < b
        @test a <= b
        @test b > a
        @test b >= a
    end

    @test CSFOrbital(1,0) < CSFOrbital(2,0) <  CSFOrbital(2,1)
    @test CSFOrbital(1,0) < CSFOrbital(2,0) <= CSFOrbital(2,0)
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
    @test nelectrons(cdef1) == 28

    cdef_combined = cdef1 * cdef2
    @test string(cdef_combined) == "1s(2,i)2s(2,i)2p(6,i)3s(2,i)3p(6,i)3d(10,i)4s(2,i)4p(6,i)4d(10,i)4f(14,i)5s(2,i)5p(6,i)5d(10,i)5f(14,i)"
    @test nelectrons(cdef_combined) == 28 + 64

    # Test `nexcitations(::CSFDefinition, ::CSFDefinition)`
    cdef3 = CSFDefinition()
    push!(cdef3, GRASP.CSFOrbital(1, 0), 2, 0)
    push!(cdef3, GRASP.CSFOrbital(2, 0), 2, 0)
    cdef4 = CSFDefinition()
    push!(cdef4, GRASP.CSFOrbital(1, 0), 2, 0)
    push!(cdef4, GRASP.CSFOrbital(2, 1), 2, 0)

    @test nexcitations(cdef3, cdef3) == 0
    @test nexcitations(cdef4, cdef4) == 0
    @test nexcitations(cdef3, cdef4) == 2

    # TODO: Test error for *?
end

import GRASP: csfdefinition, FilledOrbital, CSF
@testset "csfdefinition()" begin
    csf = let orbs = [
            FilledOrbital(1, -1, 0, 2),
            FilledOrbital(3, 2, 0, 4)
        ]
        CSF(0, true, orbs, [0, 0])
    end

    csfdef = csfdefinition(csf)
    @test nelectrons(csf) == nelectrons(csfdef)
end

end # @testset "csfs.jl"
