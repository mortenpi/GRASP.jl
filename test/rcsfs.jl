using Test
using GRASP
import AtomicLevels: @o_str, @c_str

@testset "rcsfs.jl" begin

@testset "kappa" begin
    import GRASP: kappa
    @test kappa(o"1s")  == -1
    @test kappa(o"2s")  == -1
    @test kappa(o"2p-") ==  1
    @test kappa(o"2p")  == -2
    @test kappa(o"3d-") ==  2
    @test kappa(o"3d")  == -3
end

import GRASP: CSF, Symmetries
@testset "CSF" begin
    orb1 = o"1s"  # 1s(1/2)
    orb1b = o"1s" # 1s(1/2)
    orb2 = o"3d-" # 3d(3/2)
    orb3 = o"4d"  # 4d(5/2)

    angmom = Symmetries.AngularMomentum[0, 0]

    # CSF(total2J, parity, orbs, coupled2Js)
    csf1 = CSF([orb1, orb2], [2, 4], angmom, angmom, Symmetries.even)
    @test string(csf1) == "1s(2|0|0) 3d⁻(4|0|0) | 0+"
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

import GRASP: kappa2rso
@testset "kappa2rso" begin
    @test kappa2rso(-1) == "s"
    @test kappa2rso(-2) == "p"
    @test kappa2rso(-3) == "d"
    @test kappa2rso(1)  == "p-"
    @test kappa2rso(2)  == "d-"
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
            #@test string(csf) == "1s(1|1/2|1/2) 2s(1|1/2|1) 2p(3|3/2|1/2) | 1/2-"
            @test csf.terms[1] == 1//2
            @test csf.terms[2] == 1
            @test csf.terms[3] == 1//2
        end
        let csf = csfbs[1].csfs[5]
            #@test string(csf) == "1s(1|1/2|1/2) 2s(1|1/2|1) 2p⁻(2|0|1) 2p(1|3/2|1/2) | 1/2-"
            @test csf.terms[1] == 1//2
            @test csf.terms[2] == 1
            @test csf.terms[3] == 1
            @test csf.terms[4] == 1//2
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
