using PyPlot
using LinearAlgebra
using Roots # Install by typing `] add Roots`

function euler_step(f,y0,t)
    return y0 + f(y0)*t
end

function simpson_step(f,y0,t)
    f1 = t*f(y0)
    f2 = t*f(y0 + f1/2)
    f3 = t*f(y0 + f1)
    return y0 + f1/6 + f2*4/6 + f3/6
end

function implicit_simpson_step(f,y0,t)
    f1 = t*f(y0)
    f2 = t*f(y0 + f1/2)
    return find_zero(
        y -> y0 + f1/6 + f2*4/6 + f(y)*t/6 - y,
        euler_step(f,y0,t)
    )
end

function integrate(f,y0,T,n,step)
    y = Vector{typeof(y0)}(undef,n)
    y[1] = y0
    for i = 2:n
        y[i] = step(f,y[i-1],T/(n-1))
    end
    return y
end

function convergence()
    f = y->y^2
    y0 = 1.0
    T = 0.5
    y = t-> y0/(1-y0*t)

    clf()
    n = round.(Int, 10.0.^LinRange(1,3,30))
    for (name,step) in (
        # ("explicit", simpson_step),
        ("implict", implicit_simpson_step),
    )
        error = [begin
            ỹ = integrate(f,y0,T,n, step)
            abs(y(T) - ỹ[end])
        end for n in n]
        loglog(n, error, label=name)
    end
    loglog(n, inv.(n).^2, "k--")
    legend(loc="best")
    xlabel(L"n")
    ylabel(L"|\tilde y(T) - y(T)|")
    display(gcf())
end
