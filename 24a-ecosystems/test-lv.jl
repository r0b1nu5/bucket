using PyPlot

include("lv.jl")
include("tools.jl")

n = 30

# #=
r = rand(n)
θ = rand(n)

B = rand(n,n)
#B = 2*rand(n,n) .- 1
#B = -rand(n,n) .+ 0.1
#B = gen_wER(n,.4)

# =#

x = lv(.1*rand(n),r,θ,B,.01,true,1e-6,10000)
s = round.(x[:,end];digits=2)



