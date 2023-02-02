using PyPlot

# Generates a wheel graph with 'n' nodes (including the center which is node 1) where the spokes are removed with proba p2, the tires are removed with proba p3 and triangles are added between spokes with proba p1)

function gen_rand_hyperwheel(n::Int64, p1::Float64, p2::Float64=0., p3::Float64=0., plot::Bool=false)
    A2 = zeros(n,n)
    A2[1,2] = (rand() > p2)
    A2[2,1] = A2[1,2]
    A2[2,n] = (rand() > p3)
    A2[n,2] = A2[2,n]
    for i in 3:n
        A2[1,i] = (rand() > p2)
        A2[i,1] = A2[1,i]
        A2[i,i-1] = (rand() > p3)
        A2[i-1,i] = A2[i,i-1]
    end

    e3 = zeros(3,3,3)
    e3[1,2,3] = e3[1,3,2] = e3[2,1,3] = e3[2,3,1] = e3[3,1,2] = e3[3,2,1] = 1.

    A3 = zeros(n,n,n)
    A3[[1,2,n],[1,2,n],[1,2,n]] = e3*(rand() < p1)
    for i in 3:n
        A3[[1,i-1,i],[1,i-1,i],[1,i-1,i]] = e3*(rand() < p1)
    end

    if plot
        plot_hyperwheel(A2,A3)
    end

    return A2, A3
end

function plot_hyperwheel(A2,A3)
    n = size(A2)[1]
    x = [0.; cos.((0:n-2)*2π/(n-1))]
    y = [0.; sin.((0:n-2)*2π/(n-1))]

    if A3[1,2,n] > .1
        PyPlot.fill(x[[1,2,n]], y[[1,2,n]], color="C0", alpha=.5)
    end
    for i in 3:n
        if A3[1,i-1,i] > .1
            PyPlot.fill(x[[1,i-1,i]],y[[1,i-1,i]], color="C0", alpha=.5)
        end
    end

    if A2[1,2] > .1
        PyPlot.plot(x[[1,2]], y[[1,2]], color="C0")
    end
    if A2[2,n] > .1
        PyPlot.plot(x[[2,n]], y[[2,n]], color="C0")
    end
    for i in 3:n
        if A2[1,i] > .1
            PyPlot.plot(x[[1,i]], y[[1,i]], color="C0")
        end
        if A2[i-1,i] > .1
            PyPlot.plot(x[[i-1,i]], y[[i-1,i]], color="C0")
        end
    end

    PyPlot.plot(x,y,color="C0","o")
end



