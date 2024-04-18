using PyPlot

include("lv.jl")

r = [2/3,-1.]
θ = [0.,0.]
B = [0. -4/3;1. 0.]

xs = lv(rand(2),r,θ,B,.01,true)

PyPlot.plot(xs[1,:],xs[2,:])



