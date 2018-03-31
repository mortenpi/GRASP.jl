using Base.Test
using GRASP

@testset "rcsfs.jl" begin

import GRASP: FilledOrbital, CSF
@testset "CSF" begin
    # FilledOrbital(n, kappa, total2J, nelectrons)
    orb1 = FilledOrbital(1, -1, 2, 0)
    @test string(orb1) == "1s(2~0)"

    orb2 = FilledOrbital(3, 2, 4, 0)
    @test string(orb2) == "3d-(4~0)"

    orb3 = FilledOrbital(4, -3, 4, 0)
    @test string(orb3) == "4d(4~0)"

    # CSF(total2J, parity, orbs, coupled2Js)
    csf1 = CSF(0, true, [orb1, orb2], [0, 0])
    @test string(csf1) == "1s(2~0) 3d-(4~0) ~ 0+"
    @test nelectrons(csf1) == 6

    csf2 = CSF(0, true, [orb1, orb3], [0, 0])
    @test string(csf2) == "1s(2~0) 4d(4~0) ~ 0+"
    @test nelectrons(csf2) == 6

    @test nexcitations(csf1, csf1) == 0
    @test nexcitations(csf2, csf2) == 0
    @test nexcitations(csf1, csf2) == 4
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
    let csfs = GRASP.parse_rcsf("grasp/csls/example1.c")
        @test isa(csfs, Vector{GRASP.CSFBlock})
        @test length(csfs) == 4

        @test parity(csfs[1]) == Symmetries.even
        @test angularmomentum(csfs[1]) == AngularMomentum(0)

        @test angularmomentum(csfs[2]) == AngularMomentum(1)
        @test angularmomentum(csfs[3]) == AngularMomentum(2)
        @test angularmomentum(csfs[4]) == AngularMomentum(3)
    end

    let csfs = GRASP.parse_rcsf("grasp/csls/example2.c")
        @test isa(csfs, Vector{GRASP.CSFBlock})
        @test length(csfs) == 4

        @test parity(csfs[1]) == Symmetries.odd
        @test angularmomentum(csfs[1]) == AngularMomentum(1//2)

        @test angularmomentum(csfs[2]) == AngularMomentum(3//2)
        @test angularmomentum(csfs[3]) == AngularMomentum(5//2)
        @test angularmomentum(csfs[4]) == AngularMomentum(7//2)
    end
end

end
