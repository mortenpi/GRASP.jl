using Test
using GRASP
import GRASP.Configurations: CSFOrbital, CSFDefinition, CSFDefinitionList
import AtomicLevels: @ro_str, @c_str

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

    @test maxelectrons(CSFOrbital(1,0)) == 2
    @test maxelectrons(CSFOrbital(2,0)) == 2
    @test maxelectrons(CSFOrbital(5,1)) == 6
    @test maxelectrons(CSFOrbital(5,2)) == 10
    @test maxelectrons(CSFOrbital(5,3)) == 14
    @test maxelectrons(CSFOrbital(5,4)) == 18
end

@testset "CSFDefinition" begin
    cdef1 = let cdef = CSFDefinition()
        for n=1:3, l=0:(n-1)
             push!(cdef, CSFOrbital(n, l), 2*(2l+1), 0)
        end
        cdef
    end

    cdef2 = let cdef = CSFDefinition()
        for n=4:5, l=0:min(3, n-1)
             push!(cdef, CSFOrbital(n, l), 2*(2l+1), 0)
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
    push!(cdef3, CSFOrbital(1, 0), 2, 0)
    push!(cdef3, CSFOrbital(2, 0), 2, 0)
    cdef4 = CSFDefinition()
    push!(cdef4, CSFOrbital(1, 0), 2, 0)
    push!(cdef4, CSFOrbital(2, 1), 2, 0)

    @test nexcitations(cdef3, cdef3) == 0
    @test nexcitations(cdef4, cdef4) == 0
    @test nexcitations(cdef3, cdef4) == 2

    # TODO: Test error for *?
end

import GRASP.Configurations: csfdefinition
import GRASP: CSF, Symmetries, AngularMomentum
@testset "csfdefinition()" begin
    csf = CSF([ro"1s", ro"3d-"], [2, 4], AngularMomentum[0, 0], AngularMomentum[0, 0], Symmetries.even)

    csfdef = csfdefinition(csf)
    @test nelectrons(csf) == nelectrons(csfdef)
end

@testset "parse(::CSFDefinition)" begin
    csfdef = parse(CSFDefinition, "2s(2,2) 3p(3)")
    @test isa(csfdef, CSFDefinition)
    @test length(csfdef.orbitals) == 2
    @test csfdef.orbitals[1] == CSFOrbital(2, 0)
    @test csfdef.orbitals[2] == CSFOrbital(3, 1)
    @test csfdef.nelectrons == [2, 3]
    @test csfdef.nexcitations == [2, 0]

    @test_throws ArgumentError parse(CSFDefinition, "2s(2)3p() 4p()")
    @test_throws ArgumentError parse(CSFDefinition, "2s(2)3p() 4pÖ()")
    @test_throws ArgumentError parse(CSFDefinition, "2s(2)3p() 4p a()")
    @test_throws ArgumentError parse(CSFDefinition, "2s 3pÖ()()")
    @test_throws ArgumentError parse(CSFDefinition, "2s 3pÖÖ(")
    @test_throws ArgumentError parse(CSFDefinition, "2s 3pÖÖ(())")
end

import GRASP.Configurations: parse_l
@testset "parse_l" begin
    @test parse_l("s")   == 0
    @test parse_l("s  ") == 0
    @test parse_l("p")   == 1
    @test parse_l("h")   == 5
end

import GRASP.Configurations: parse_j
@testset "parse_j" begin
    @test parse_j("s")  == -1
    @test parse_j("p-") ==  1
    @test parse_j("p")  == -2
    @test parse_j("d-") ==  2
    @test parse_j("d")  == -3
    @test parse_j("f-") ==  3
end

end # @testset "csfs.jl"
