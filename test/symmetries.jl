module SymmetriesTests
using Base.Test
using GRASP.Symmetries

#
# Parity
# ------------------------------------------------------------------------------------------
@test isdefined(Symmetries, :even)
@test isdefined(Symmetries, :odd)

@test Parity('+') == Symmetries.even
@test Parity('-') == Symmetries.odd
@test_throws ArgumentError Parity(' ')

@test Parity("+") == Symmetries.even
@test Parity("-") == Symmetries.odd
@test_throws ArgumentError Parity("")
@test_throws ArgumentError Parity("+-")
@test_throws ArgumentError Parity("a")
@test_throws ArgumentError Parity("α")

@test string(Symmetries.even) == "+"
@test string(Symmetries.odd) == "-"

#
# Angular momentum
# ------------------------------------------------------------------------------------------
@test_throws ArgumentError AngularMomentum(-2)

@test string(AngularMomentum(0))  == "0"
@test string(AngularMomentum(1))  == "1/2"
@test string(AngularMomentum(2))  == "1"
@test string(AngularMomentum(3))  == "3/2"
@test string(AngularMomentum(10)) == "5"
@test string(AngularMomentum(21)) == "21/2"

#
# Angular symmetries
# ------------------------------------------------------------------------------------------
@test string(AngularSymmetry(AngularMomentum(0), Symmetries.even)) == "0+"
@test string(AngularSymmetry(AngularMomentum(1), Symmetries.odd))  == "1/2-"

@test string(AngularSymmetry(0, Symmetries.even))  == "0+"
@test string(AngularSymmetry(1, Symmetries.odd))   == "1/2-"

@test string(AngularSymmetry(AngularMomentum(0), '-'))  == "0-"
@test string(AngularSymmetry(AngularMomentum(1), "+"))  == "1/2+"

@test string(AngularSymmetry(0, true))  == "0+"
@test string(AngularSymmetry(1, false)) == "1/2-"
@test string(AngularSymmetry(2, '-'))   == "1-"
@test string(AngularSymmetry(3, "+"))   == "3/2+"

end
