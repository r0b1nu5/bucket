using Statistics

include("kuramoto.jl")
include("plot_circ.jl")

n = 20
P0 = 5.
maxiter = 100000
eps = 1e-7
dt = 1e-3

A = 1 .- diagm(0 => ones(n))
D = (n-1)*diagm(0 => ones(n))
L = D - A

P = P0*(2*rand(n) .- 1)
P .-= mean(P)
th0 = 2pi*rand(n)

ths,dth,it = kuramoto(L,P,th0,true,true,maxiter,eps,dt)
plot_circ_log(ths)
plot_circ_lin(ths)


