using DelimitedFiles, LinearAlgebra

include("L2B.jl")

n = 10
m = 24
m2 = 12

L = readdlm("ntw_data/ntw10_L.csv",',')
include("ntw_data/ntw10_cycles.jl")

b,w = L2B(L)

B1 = b.*(b .> 0)
B2 = -b.*(b .< 0)

Bout = [B1 B2]
Bin = [B2 B1]

B = Bout - Bin

ω = zeros(n)
γ = π/4

function h(x::Float64, α::Float64=.1)
	return sin(x - α) + sin(α)
end

function hi(f::Float64, α::Float64=.1)
	return asin(f - sin(α)) + α
end

function H(f::Float64, α::Float64=.1)
	return h(-hi(f,α),α)
end

I2 = diagm(0 => ones(m2))
L2 = .5*diagm(0 => ones(m2))

P = I2 - L2*b'*pinv(b*L2*b')*b

function Tu(f::Array{Float64,1}, u::Array{Int64,1}, P::Array{Float64,2}, C::Array{Float64,2})
	gif = [hi(f[i]) for i in 1:m2]

	x = gif - P*L2*(gif - 2*π*pinv(C)*u)

	return [[h(x[i]) for i in 1:m2];[h(-x[i]) for i in 1:m2]]
end







