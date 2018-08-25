using Test
using GRASP.Symmetries

@testset "Symmetries" begin

#
# Parity
# ------------------------------------------------------------------------------------------
@test isdefined(Symmetries, :even)
@test isdefined(Symmetries, :odd)

@test Parity('+') == Symmetries.even
@test Parity('-') == Symmetries.odd
@test_throws ArgumentError Parity(' ')

@test string(Symmetries.even) == "+"
@test string(Symmetries.odd)  == "-"

@test parse(Parity, "+")   == Symmetries.even
@test parse(Parity, "-")   == Symmetries.odd
@test parse(Parity, " + ") == Symmetries.even
@test parse(Parity, "\t+") == Symmetries.even

@test_throws Meta.ParseError parse(Parity, "")
@test_throws Meta.ParseError parse(Parity, " ")
@test_throws Meta.ParseError parse(Parity, "++")
@test_throws Meta.ParseError parse(Parity, "α")

#
# Angular momentum
# ------------------------------------------------------------------------------------------
@test_throws ArgumentError AngularMomentum(-2)
@test_throws ArgumentError AngularMomentum(-1//2)

@test string(AngularMomentum(0))     == "0"
@test string(AngularMomentum(1//2))  == "1/2"
@test string(AngularMomentum(2//2))  == "1"
@test string(AngularMomentum(3//2))  == "3/2"
@test string(AngularMomentum(5))     == "5"
@test string(AngularMomentum(21//2)) == "21/2"

@test parse(AngularMomentum, "")     == AngularMomentum(0)
@test parse(AngularMomentum, "  ")   == AngularMomentum(0)
@test parse(AngularMomentum, "0")    == AngularMomentum(0)
@test parse(AngularMomentum, "1")    == AngularMomentum(1)
@test parse(AngularMomentum, "22")   == AngularMomentum(22)
@test parse(AngularMomentum, "1/2")  == AngularMomentum(1//2)
@test parse(AngularMomentum, "11/2") == AngularMomentum(11//2)
@test parse(AngularMomentum, "10/2") == AngularMomentum(10//2)


#
# Angular symmetries
# ------------------------------------------------------------------------------------------
@test string(AngularSymmetry(AngularMomentum(0), Symmetries.even))   == "0+"
@test string(AngularSymmetry(AngularMomentum(1//2), Symmetries.odd)) == "1/2-"

@test string(AngularSymmetry(0, Symmetries.even))   == "0+"
@test string(AngularSymmetry(1//2, Symmetries.odd)) == "1/2-"

@test string(AngularSymmetry(AngularMomentum(0), '-'))    == "0-"
@test string(AngularSymmetry(AngularMomentum(1//2), '+')) == "1/2+"

@test string(AngularSymmetry(0, true))     == "0+"
@test string(AngularSymmetry(1//2, false)) == "1/2-"
@test string(AngularSymmetry(1, '-'))      == "1-"
@test string(AngularSymmetry(3//2, '+'))   == "3/2+"

end # @testset "Symmetries"
