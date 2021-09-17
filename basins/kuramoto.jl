using PyPlot, DelimitedFiles, LinearAlgebra, SparseArrays

include("L2B.jl")

# Runs RK4 simulation of the Kuramoto model on a cycle.
# Returns the final state.
function kuramoto(L::SparseMatrixCSC{Float64, Int64}, ω::Array{Float64,1}, θ0::Array{Float64,1}, verb::Bool=false, h::Float64=.1, max_iter::Int64=10000, tol::Float64=1e-6)
	B,w,Bt = L2B(L)
	W = spdiagm(0 => w)
	BW = B*W

	n,m = size(B)

	θ = θ0
	iter = 0
	err = 1000.

	while err > tol && iter < max_iter
		iter += 1
		if verb && iter%1000 == 0
			@info "iter = $iter"
		end

		Btθ = Bt*θ
		k1 = ω - BW*sin.(Btθ)
		k2 = ω - BW*sin.(Btθ + h/2*Bt*k1)
		k3 = ω - BW*sin.(Btθ + h/2*Bt*k2)
		k4 = ω - BW*sin.(Btθ + h*Bt*k3)

		dθ = (k1 + 2*k2 + 2*k3 + k4)/6

		θ += h*dθ

		θ = mod.(θ .+ π,2π) .- π

		err = maximum(abs.(dθ))
	end

	return θ,iter
end

# Runs RK4 simulation of the Kuramoto model on a cycle. 
# Returns the whole time series.
function kuramoto_series(L::SparseMatrixCSC{Float64, Int64}, ω::Array{Float64,1}, θ0::Array{Float64,1}, h::Float64=.1, max_iter::Int64=10000, tol::Float64=1e-6)
	B,w,Bt = L2B(L)
	W = spdiagm(0 => w)
	BW = B*W

	n,m = size(B)

	θ = θ0
	θs = Array{Float64,2}(undef,n,0)
	iter = 0
	err = 1000.

	while err > tol && iter < max_iter
		iter += 1
		if iter%1000 == 0
			@info "iter = $iter"
		end

		θs = [θs θ]

		Btθ = Bt*θ
		k1 = ω - BW*sin.(Btθ)
		k2 = ω - BW*sin.(Btθ + h/2*Bt*k1)
		k3 = ω - BW*sin.(Btθ + h/2*Bt*k2)
		k4 = ω - BW*sin.(Btθ + h*Bt*k3)

		dθ = (k1 + 2*k2 + 2*k3 + k4)/6

		θ += h*dθ
		θ = mod.(θ .+ π,2π) .- π

		err = maximum(abs.(dθ))
	end

	return θs
end

# Computes the winding number of θ around cycle cyc.
function winding(θ::Array{Float64,1}, cyc::Array{Int64,1})
	dθ = mod.(θ[[cyc[2:end];cyc[1]]] - θ[cyc] .+ π,2π) .- π

	return round(Int64,sum(dθ)/(2π))
end

# Computes the winding vector around cycles in cyc.
function winding(θ::Array{Float64,1}, cyc::Array{Array{Int64,1},1})
	return [winding(θ,c) for c in cyc]
end

# Computes the distance between θ1 and θ2 on the torus.
function dist(θ1::Array{Float64,1}, θ2::Array{Float64,1})
	return torus_norm(θ1 - θ2)
end

# Computes the torus norm an angle vector, modulo constant shift. 
# Makes sure that the vector has components between -π and π, and is orhtogonal to the constant vector.
function torus_norm(θ::Array{Float64,1}, max_iter::Int64=100)
	θ1 = θ
	θ2 = θ
	no = 1000.
	iter = 0
	while no > 1e-8 && iter < max_iter
		iter += 1
		θ1 = θ2
		θt = mod.(θ1 .+ π,2π) .- π
		θ2 = θt .- mean(θt)
		no = norm(θ1 - θ2)
	end

	if iter == max_iter
		@info "torus_norm: Did not converge..."
	end

	return norm(θ2)
end

# Computes the distance between an angle vector θ and the splay state with winding number q.
function dist2splay(θ::Array{Float64,1}, q::Int64)
	n = length(θ)

	sq = Array(2π*q*(0:n-1)/n)

	no = Inf

	for i in 1:n
		no = min(no,torus_norm(θ - sq[[(i:n);(1:i-1)]]))
	end

	return no
end
	
# Computes the winding number of the basin where θ lies.
function basin_test(θ::Vector{Float64}, L::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, ω::Vector{Float64})
	θf = kuramoto(L,ω,θ)
	
	return winding(θf)
end




