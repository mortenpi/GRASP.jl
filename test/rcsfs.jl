using Test
using GRASP

@testset "rcsfs.jl" begin

import GRASP: RelativisticOrbital
@testset "RelativisticOrbital" begin
    orbs = [
        RelativisticOrbital(1, -1), # 1s  / 1s(1/2)

        RelativisticOrbital(2, -1), # 2s  / 2s(1/2)
        RelativisticOrbital(2,  1), # 2p- / 2p(1/2)
        RelativisticOrbital(2, -2), # 2p  / 2p(3/2)

        RelativisticOrbital(3, -1), # 3s  / 3s(1/2)
        RelativisticOrbital(3,  1), # 3p- / 3p(1/2)
        RelativisticOrbital(3, -2), # 3p  / 3p(3/2)
        RelativisticOrbital(3,  2), # 3d- / 3d(3/2)
        RelativisticOrbital(3, -3), # 3d  / 3d(5/2)
    ]

    @test string(orbs[1]) == "1s"

    @test string(orbs[2]) == "2s"
    @test string(orbs[3]) == "2p-"
    @test string(orbs[4]) == "2p"

    @test string(orbs[5]) == "3s"
    @test string(orbs[6]) == "3p-"
    @test string(orbs[7]) == "3p"
    @test string(orbs[8]) == "3d-"
    @test string(orbs[9]) == "3d"

    for orb in orbs
        @test orb == orb
    end

    for i=1:length(orbs), j=(i+1):length(orbs)
        a, b = orbs[i], orbs[j]
        @test a != b
        @test a < b
        @test a <= b
        @test b > a
        @test b >= a
    end

    @test orbs[1] < orbs[2] < orbs[3] < orbs[4]
    @test orbs[1] <= orbs[2] <= orbs[3] <= orbs[4]

    # angularmomentum(::RelativisticOrbital)
    @test angularmomentum(orbs[1]) == GRASP.AngularMomentum(1//2)
    @test angularmomentum(orbs[2]) == GRASP.AngularMomentum(1//2)
    @test angularmomentum(orbs[3]) == GRASP.AngularMomentum(1//2)
    @test angularmomentum(orbs[4]) == GRASP.AngularMomentum(3//2)
    @test angularmomentum(orbs[5]) == GRASP.AngularMomentum(1//2)
    @test angularmomentum(orbs[6]) == GRASP.AngularMomentum(1//2)
    @test angularmomentum(orbs[7]) == GRASP.AngularMomentum(3//2)
    @test angularmomentum(orbs[8]) == GRASP.AngularMomentum(3//2)
    @test angularmomentum(orbs[9]) == GRASP.AngularMomentum(5//2)

    @test maxelectrons(orbs[1]) == 2

    @test maxelectrons(orbs[2]) == 2
    @test maxelectrons(orbs[3]) == 2
    @test maxelectrons(orbs[4]) == 4

    @test maxelectrons(orbs[5]) == 2
    @test maxelectrons(orbs[6]) == 2
    @test maxelectrons(orbs[7]) == 4
    @test maxelectrons(orbs[8]) == 4
    @test maxelectrons(orbs[9]) == 6
end

import GRASP: CSF, Symmetries
@testset "CSF" begin
    orb1 = RelativisticOrbital(1, -1) # 1s  / 1s(1/2)
    orb1b = RelativisticOrbital(1, -1) # 1s  / 1s(1/2)
    orb2 = RelativisticOrbital(3,  2) # 3d- / 3d(3/2)
    orb3 = RelativisticOrbital(4, -3) # 4d  / 4d(5/2)

    angmom = Symmetries.AngularMomentum[0, 0]

    # CSF(total2J, parity, orbs, coupled2Js)
    csf1 = CSF([orb1, orb2], [2, 4], angmom, angmom, Symmetries.even)
    @test string(csf1) == "1s(2|0|0) 3d-(4|0|0) | 0+"
    @test nelectrons(csf1) == 6

    csf1b = CSF([orb1b, orb2], [2, 4], angmom, angmom, Symmetries.even)
    @test csf1 == csf1b

    csf2 = CSF([orb1, orb3], [2, 4], angmom, angmom, Symmetries.even)
    @test string(csf2) == "1s(2|0|0) 4d(4|0|0) | 0+"
    @test nelectrons(csf2) == 6

    @test nexcitations(csf1, csf1) == 0
    @test nexcitations(csf2, csf2) == 0
    @test nexcitations(csf1, csf2) == 4

    # Test the ::CSF iterator
    @test length(csf1) == 2
    let cs = collect(csf1), zero_angmom = Symmetries.AngularMomentum(0)
        @test cs[1][1] == orb1
        @test cs[1][2] == 2
        @test cs[1][3] == zero_angmom
        @test cs[1][4] == zero_angmom

        @test cs[2][1] == orb2
        @test cs[2][2] == 4
        @test cs[2][3] == zero_angmom
        @test cs[2][4] == zero_angmom
    end
end

import GRASP: parse_l
@testset "parse_l" begin
    @test parse_l("s")   == 0
    @test parse_l("s  ") == 0
    @test parse_l("p")   == 1
    @test parse_l("h")   == 5
end

import GRASP: parse_j
@testset "parse_j" begin
    @test parse_j("s")  == -1
    @test parse_j("p-") ==  1
    @test parse_j("p")  == -2
    @test parse_j("d-") ==  2
    @test parse_j("d")  == -3
    @test parse_j("f-") ==  3
end

import GRASP: kappa2rso
@testset "kappa2rso" begin
    @test kappa2rso(-1) == "s"
    @test kappa2rso(-2) == "p"
    @test kappa2rso(-3) == "d"
    @test kappa2rso(1)  == "p-"
    @test kappa2rso(2)  == "d-"
end

import GRASP: parse_orbital
@testset "parse_orbital" begin
    @test parse_orbital("1s")  == (1, -1)
    @test parse_orbital("2p-") ==  (2, 1)
    @test parse_orbital("33p") == (33, -2)
    @test parse_orbital("142d-") ==  (142, 2)
    @test parse_orbital(" 124214  d") == (124214, -3) # Maybe should throw an error?
    @test parse_orbital("   1f-") ==  (1, 3)
end

import GRASP: Symmetries, AngularMomentum, angularmomentum, parity
@testset "parse_rcsf" begin
    let csfbs = GRASP.parse_rcsf(joinpath(@__DIR__, "grasp/csls/example1.c"))
        @test isa(csfbs, Vector{GRASP.CSFBlock})
        @test length(csfbs) == 4

        @test parity(csfbs[1]) == Symmetries.even
        @test angularmomentum(csfbs[1]) == AngularMomentum(0)

        @test angularmomentum(csfbs[2]) == AngularMomentum(1)
        @test angularmomentum(csfbs[3]) == AngularMomentum(2)
        @test angularmomentum(csfbs[4]) == AngularMomentum(3)
    end

    let csfbs = GRASP.parse_rcsf(joinpath(@__DIR__, "grasp/csls/example2.c"))
        @test isa(csfbs, Vector{GRASP.CSFBlock})
        @test length(csfbs) == 4

        @test parity(csfbs[1]) == Symmetries.odd
        @test angularmomentum(csfbs[1]) == AngularMomentum(1//2)

        @test angularmomentum(csfbs[2]) == AngularMomentum(3//2)
        @test angularmomentum(csfbs[3]) == AngularMomentum(5//2)
        @test angularmomentum(csfbs[4]) == AngularMomentum(7//2)

        # Test that we interpret the in-shell and inter-shell couplings correctly
        let csf = csfbs[1].csfs[1]
            @test string(csf) == "1s(1|1/2|1/2) 2s(1|1/2|1) 2p(3|3/2|1/2) | 1/2-"
            @test csf.csfcouplings[1] == AngularMomentum(1//2)
            @test csf.csfcouplings[2] == AngularMomentum(1)
            @test csf.csfcouplings[3] == AngularMomentum(1//2)
        end
        let csf = csfbs[1].csfs[5]
            @test string(csf) == "1s(1|1/2|1/2) 2s(1|1/2|1) 2p-(2|0|1) 2p(1|3/2|1/2) | 1/2-"
            @test csf.csfcouplings[1] == AngularMomentum(1//2)
            @test csf.csfcouplings[2] == AngularMomentum(1)
            @test csf.csfcouplings[3] == AngularMomentum(1)
            @test csf.csfcouplings[4] == AngularMomentum(1//2)
        end
    end

    let csfbs = GRASP.parse_rcsf(joinpath(@__DIR__, "grasp/csls/example3.c"))
        @test isa(csfbs, Vector{GRASP.CSFBlock})
        @test length(csfbs) == 6

        @test parity(csfbs[1]) == Symmetries.even
        @test angularmomentum(csfbs[1]) == AngularMomentum(1//2)
        @test angularmomentum(csfbs[2]) == AngularMomentum(3//2)
        @test angularmomentum(csfbs[3]) == AngularMomentum(5//2)
        @test angularmomentum(csfbs[4]) == AngularMomentum(7//2)
        @test angularmomentum(csfbs[5]) == AngularMomentum(9//2)
        @test angularmomentum(csfbs[6]) == AngularMomentum(11//2)
    end
end

end # @testset "rcsfs.jl"
