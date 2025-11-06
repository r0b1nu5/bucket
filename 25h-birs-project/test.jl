using PyPlot

include("sivo.jl")

n = 10
p = .6

s0 = .99*ones(n)
x0 = .01*ones(n)
v0 = 0*ones(n)
o0 = rand(n) .- .5

A = rand(n,n)
A = (A + A').*(1 .- diagm(0 => ones(n))) .- 1

B = (rand(n,n) .< p).*(1 .- diagm(0 => ones(n)))

γ = .1

ρ = .3

S,X,V,O = sivo(s0,x0,v0,o0,A,B,γ,ρ,.1)

for i in 1:n
    figure("fig")
    subplot(2,2,1)
    PyPlot.plot(S[i,:])
    ylabel("s")
    subplot(2,2,2)
    PyPlot.plot(X[i,:])
    ylabel("x")
    subplot(2,2,3)
    PyPlot.plot(V[i,:])
    ylabel("v")
    subplot(2,2,4)
    PyPlot.plot(O[i,:])
    ylabel("o")
end





