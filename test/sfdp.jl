using Base.Test
using JLab

# #############################################################################
#Test significant figures
println("########## Testing \"to_sf\" ##########")

#Non-positive sigfig error
println("Testing for error on non-positive argument num_sf ...")
@test_throws to_sf(1,0)
@test_throws to_sf(1,-1)

#Input types
println("Testing input types ...")
@test to_sf(1,3) == "1.00"
@test to_sf(1.0,3) == "1.00"
@test to_sf(1e0,3) == "1.00"
@test to_sf(pi,3) == "3.14"

#Zero
println("Testing 0 as input ...")
@test to_sf(0,3) == "0"

#Rounding
println("Testing rounding ...")
@test to_sf(3.141593,6) == "3.14159"
@test to_sf(3.141593,5) == "3.1416"
@test to_sf(3.141593,4) == "3.142"
@test to_sf(3.141593,3) == "3.14"
@test to_sf(3.141593,2) == "3.1"
@test to_sf(3.141593,1) == "3"

#Exponents
println("Testing exponent notation ...")
@test to_sf(3.142e4,5) == "3.1420e4"
@test to_sf(3.142e4,4) == "3.142e4"

@test to_sf(3.142e3,5) == "3.1420e3"
@test to_sf(3.142e3,4) == "3142"
@test to_sf(3.142e3,3) == "3140"
@test to_sf(3.142e3,2) == "3100"
@test to_sf(3.142e3,1) == "3000"
@test to_sf(3.142e3,4,force=true) == "3.142e3"
@test to_sf(3.142e3,3,force=true) == "3.14e3"
@test to_sf(3.142e3,2,force=true) == "3.1e3"
@test to_sf(3.142e3,1,force=true) == "3e3"

@test to_sf(3.142e2,4) == "3.142e2"
@test to_sf(3.142e2,3) == "314"
@test to_sf(3.142e2,2) == "310"
@test to_sf(3.142e2,1) == "300"
@test to_sf(3.142e2,3,force=true) == "3.14e2"
@test to_sf(3.142e2,2,force=true) == "3.1e2"
@test to_sf(3.142e2,1,force=true) == "3e2"

@test to_sf(3.142e1,3) == "3.14e1"
@test to_sf(3.142e1,2) == "31"
@test to_sf(3.142e1,1) == "30"
@test to_sf(3.142e1,2,force=true) == "3.1e1"
@test to_sf(3.142e1,1,force=true) == "3e1"

@test to_sf(3.142,2) == "3.1"
@test to_sf(3.142,1,force=true) == "3"
@test to_sf(3.142,1) == "3"
@test to_sf(3.142,1,force=true) == "3"

@test to_sf(3.142e-1,3) == "0.314"
@test to_sf(3.142e-1,3,force=true) == "3.14e-1"

@test to_sf(3.142e-2,3) == "3.14e-2"

#Exponent format
println("Testing exponent format ...")
@test to_sf(3.14e-2,3,expformat="\\times10^{%}") == "3.14\\times10^{-2}"

println("Function \"to_sf\" passes tests.")

# #############################################################################
#Test decimal places
println("########## Testing \"to_dp\" ##########")

#Input types
println("Testing input types ...")
@test to_dp(1,2) == "1.00"
@test to_dp(1.0,2) == "1.00"
@test to_dp(1e0,2) == "1.00"
@test to_dp(pi,2) == "3.14"

#Rounding
println("Testing rounding ...")
@test to_dp(3.141593,7) == "3.1415930"
@test to_dp(3.141593,6) == "3.141593"
@test to_dp(3.141593,5) == "3.14159"
@test to_dp(3.141593,4) == "3.1416"
@test to_dp(3.141593,3) == "3.142"

#Decimal places
println("Testing decimal places ...")
@test to_dp(3.142e4,2) == "31420.00"
@test to_dp(3.142e3,2) == "3142.00"
@test to_dp(3.142e3,0) == "3142"
@test to_dp(3.142e3,-1) == "3140"
@test to_dp(3.142e3,-3) == "3000"
@test to_dp(3.142e3,-4) == "0"

println("Function \"to_dp\" passes tests.")
