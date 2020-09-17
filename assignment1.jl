using PyPlot
using LinearAlgebra
using Random

#############
# Question 3

function expm1_errors()
    # TODO: Sample x over an appropriate range
    # Hint: Have a look at Julia's `LinRange` function.
    n = 100000
    x = LinRange(-1,1,n)
    f = x -> exp(x) - 1
    relerror = @. abs(f(x) - f(big(x))) / abs(f(big(x)))

    clf()
    # TODO: Plot `(x,relerror)` to show `relerror = O(inv(x))`.
    # plot(x, relerror)
    plot((1 ./ x)[1:Int(n/2)], relerror[1:Int(n/2)])
    plot((1 ./ x)[Int(n/2)+1:end], relerror[Int(n/2)+1:end])
    display(gcf())

    println()
    println(relerror)
    println("P(relerror < 1e-4) = ", sum(relerror .< 1e-4)/n)
    println("P(relerror > 1e-1) = ", sum(relerror .> 1e-1)/n)
end


#############
# Question 4

function linear_system_error()
    # Assemble the linear system
    Random.seed!(42)
    n = 20
    xx = LinRange(-1,1,n)
    A = xx.^(0:n-1)'
    b = rand(n)

    # Compute the solution with `Float64` and `BigFloat` accuracy
    x = x_F64 = A\b
    x_big = big.(A) \ big.(b)

    # Compute estimated and exact relative errors
    sigma = svdvals(A)

    C = sigma[1] * (1 / sigma[end])

    C_ = sigma[1] * (1 / sigma[end])
    C = C_ / (1-C_*eps())

    relerror = norm(x - x_big) / norm(x_big)

    println("Estimated relative error: ", round(C * eps(), sigdigits=3))
    println("    Exact relative error: ", round(relerror, sigdigits=3))
end
