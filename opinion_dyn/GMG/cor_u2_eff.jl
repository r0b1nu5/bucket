using PyPlot, LinearAlgebra

include("big_rand.jl")

n = 201

m1 = -.5
s1 = .2
m2 = .5
s2 = .2

ep = .2

x0 = big_rand(n,m1,s1,m2,s2)

A = Float64.((0 .< abs.(repeat(x0,1,n) - repeat(x0',n,1)) .< ep))
L = diagm(0 => vec(sum(A,dims=1))) - A

xf = consensus(L,x0,zeros(n),ones(n),ones(n))

d0 = 0.
for i in 1:n
	d0 += max(xf[i],0.)
	global d0,xf
end

ds = zeros(n)
for i in 1:n
	xr = zeros(n)
	xr[i] = -1.
	x = consensus(L,x0,xr,ones(n),ones(n))
	for j in 1:n
		ds[i] += max(x[j],0.)
	end
	global L,x0,n,ds
end

ei = eigen(L)
us = ei.vectors
ls = ei.values

i = 1
while abs(ls[i]) < 1e-8
	global i,us,d0,ds,ep,m1,m2
	i += 1
	PyPlot.plot(us[:,i],d0 .- ds,"x",label="i = $i, r = $(round(cor(us[:,i],d0 .- ds),digits=2)), ε = $ep, δ = $(m2-m1)")
end





