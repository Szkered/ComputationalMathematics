using PyPlot
using LinearAlgebra

function laplacian_1d(n)
    l =  (n)^2 * Tridiagonal(
        fill( 1.0,n-1), # subdiagonal
        fill(-2.0,n),   # diagonal
        fill( 1.0,n-1)  # superdiagonal
    )
    l[n,n-1] = (n)^2 * 2
    return l
end

function solve_poisson_1d(f, n)
    x = LinRange(0,1,n+1)[2:end]
    Δ = laplacian_1d(n)
    return x, -Δ\f.(x)
end

function example_1d()
    f = x -> π^2*sin(0.5*π*x)
    uref = x -> 4*sin(0.5*π*x)
    x,u = solve_poisson_1d(f,4)

    clf()
    xx = LinRange(0,1,1000)
    plot(xx,uref.(xx), "k-", label="reference")
    plot([0;x;1],[0;u;0], label="FD")
    legend(frameon=false)
    xlabel(L"x")
    ylabel(L"u(x)")
    display(gcf())
end

function convergence()
    # Define problem and solution
    f = x -> π^2 * sin(0.5*π*x)
    u = x -> 4*sin(0.5*π*x)

    # Compute errors
    n = 2 .^ (1:15)
    error = [begin
             x,ũ = solve_poisson_1d(f,n)
             norm(ũ .- u.(x), 2)/sqrt(n+1)
             end for n in n]

    # Plot
    clf()
    loglog(n, error, label=L"\|u - u_n\|_{2,n}")
    loglog(n, n.^-2, "k--", label=L"O(n^{-2})")
    xlabel(L"n")
    legend(frameon=false)
    display(gcf())

end


function conjugate_gradients(A,b,m)
    x = 0.0
    p = r = b
    for k = 1:m
        a = (transpose(r)*r) / (transpose(p)*A*p)
        x = x .+ a*p
        r_k = r - a*A*p
        β = (transpose(r_k)*r_k) / (transpose(r)*r)
        p = r_k + β*p
        r = r_k
    end
    return x
end

function test()
    n = 10
    m = 12
    A = rand(n,n)
    A = transpose(A)*A  # make A positive semidefinite
    b = rand(n)
    x = conjugate_gradients(A,b,m)

    @assert isapprox(A*x, b), "Test failed"
    println("Test passed")
end
