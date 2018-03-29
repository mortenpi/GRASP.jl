module ParitiesTests
using Base.Test
using GRASP.Parities

@test isdefined(Parities, :even)
@test isdefined(Parities, :odd)

@test Parity('+') == Parities.even
@test Parity('-') == Parities.odd
@test_throws ArgumentError Parity(' ')

@test Parity("+") == Parities.even
@test Parity("-") == Parities.odd
@test_throws ArgumentError Parity("")
@test_throws ArgumentError Parity("+-")
@test_throws ArgumentError Parity("a")
@test_throws ArgumentError Parity("α")

@test string(Parities.even) == "+"
@test string(Parities.odd) == "-"

end
