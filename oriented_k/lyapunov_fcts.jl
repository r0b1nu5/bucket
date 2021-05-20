using PyPlot, LinearAlgebra

include("tools.jl")

function V1(θs::Union{Array{Float64,1},Array{Float64,2}}, L::Array{Float64,2})
	B,S,T,w = L2B_dir(L)

	V = vec(w'*(1 .- cos.(B'*θs)))

	figure("V")
	PyPlot.plot(V,label="V1")

	return V
end


function V2(θs::Union{Array{Float64,1},Array{Float64,2}})
	n = size(θs)[1]

	B,S,T,w = L2B_dir(n*diagm(0 => ones(n)) - ones(n,n))

	V = vec(w'*(1 .- cos.(B'*θs)))

	figure("V")
	PyPlot.plot(V,label="V2")

	return V
end

