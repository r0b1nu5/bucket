using PyPlot, LinearAlgebra

include("sivo.jl")

n = 10
p = .3
α = 1.
β = .2
h = .01

s0 = .99*ones(n)
x0 = .01*ones(n)
v0 = 0*ones(n)
o0 = rand(n) .- .5

A = α*(2*rand(n,n) .- 1)
A = (A + A').*(1 .- diagm(0 => ones(n))) 
a = abs.(A)*ones(n)
AA = A - diagm(0 => a)

B = β*(rand(n,n) .< p).*(1 .- diagm(0 => ones(n)))

γ = .1

ρ = .1

S,X,V,O = sivo(s0,x0,v0,o0,A,B,γ,ρ,h,10000)

for i in 1:n
    figure("fig",(10,8))
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

Jxx = get_Jxx(S[:,end],X[:,end],V[:,end],O[:,end],A,B,γ)
λs = eigvals(Jxx)

figure("eigenvalues")
PyPlot.plot(real.(λs),imag.(λs),"o")
xlabel("Re(λ)")
ylabel("Im(λ)")



