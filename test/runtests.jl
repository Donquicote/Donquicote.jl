using Donquicote, Test, BenchmarkTools
using Roots, ForwardDiff

f(x) = exp(x) - x^4
∂f(x) = exp(x) - 3*x^4
∂2f(x) = exp(x) - 2*x^3
D(f) = x -> ForwardDiff.derivative(f,float(x))

println("Bisection (built-in):")
@btime find_zero(f, (1,3), Bisection(); xatol = 1e-7)
println("FalsePosition (build-in):")
@btime find_zero(f, (1,3), FalsePosition(), xatol = 1e-7)
println("Newton (built-in):")
@btime find_zero((f, D(f)), 2, Roots.Newton(); xatol = 1e-7)
println("Newton (custom):")
@btime newton(f, D(f), 2)
println("Halley (built-in):")
@btime find_zero((f, D(f), D(D(f))), 2, Roots.Halley(); xatol = 1e-7)
println("Halley (custom):")
@btime haley(f, D(f), D(D(f)), 2)
println("Schroeder (built-in):")
@btime find_zero((f, D(f), D(D(f))), 2, Roots.Schroeder(); xatol = 1e-7)
println("Schroeder (custom):")
@btime schroeder(f, D(f), D(D(f)), 2)