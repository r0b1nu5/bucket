using PyPlot, LinearAlgebra

include("kuramoto.jl")
include("cycle.jl")

n = 10
L = cycle(n)
ω = zeros(n)
h = .001
T = 20000

N = 1000
d0 = 1e-3

maxd = Vector{Float64}()


for ξ in 1:N
	global L,ω,d0,n

	θ0 = 2π*rand(n) .- π
	v = 2*rand(n) .- 1
	v ./= norm(v)
	θ1 = θ0 + d0*v

	θs0 = kuramoto(L,ω,θ0,false,h,T,-1.)
	θs1 = kuramoto(L,ω,θ1,false,h,T,-1.)

	δ = θs0 - θs1

	d = sum(δ.^2,dims=1)

	push!(maxd,maximum(d))
end

PyPlot.hist(maxd,50)




