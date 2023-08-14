using PyPlot

include("tools.jl")

# w: edge weights w[i] is the weight of (i,i+1).
function cycle_kuramoto(θ0::Vector{Float64}, ω::Vector{Float64}, w::Vector{Float64}, h::Float64=.01, max_iter::Int64=1000, thr::Float64=1e-6)
	θ = θ0

	t = 0
	corr = 1000.
	while t < max_iter && (corr > thr || t < 10)
		t += 1
#		if t%100 == 0
#			@info "$t/$max_iter"
#		end
		
		k1 = f_cycle_k(θ,ω,w)
		k2 = f_cycle_k(θ+k1*h/2,ω,w)
		k3 = f_cycle_k(θ+k2*h/2,ω,w)
		k4 = f_cycle_k(θ+k3*h,ω,w)

		dθ = (k1+2*k2+2*k3+k4)./6
		θ += h*dθ
		corr = norm(dθ,Inf)
	end

	return θ
end

function cycle_kuramoto(θ0::Vector{Float64}, ω::Vector{Float64}, h::Float64=.01, max_iter::Int64=1000, thr::Float64=1e-6)
	return cycle_kuramoto(θ0,ω,ones(length(θ0)),h,max_iter,thr)
end

function f_cycle_k(θ::Vector{Float64}, ω::Vector{Float64}, w::Vector{Float64})
	n = length(θ)

	return ω - w.*sin.(θ - [θ[2:end];θ[1]]) - [w[n];w[1:n-1]].*sin.(θ - [θ[n];θ[1:n-1]])
end

function f_cycle_k(θ::Vector{Float64}, ω::Vector{Float64}, h::Float64=.01, max_iter::Int64=1000, thr::Float64=1e-6)
	return f_cycle_k(θ,ω,ones(length(θ)),h,max_iter,thr)
end


	

# w: edge weights w[i] is the weight of (i,i+1).
function dir_cycle_kuramoto(θ0::Vector{Float64}, ω::Vector{Float64}, w::Vector{Float64}, h::Float64=.01, max_iter::Int64=1000, thr::Float64=1e-6)
	θ = θ0

	t = 0
	corr = 1000.
	while t < max_iter && (corr > thr || t < 10)
		t += 1
#		if t%100 == 0
#			@info "$t/$max_iter"
#		end
		
		k1 = f_dir_cycle_k(θ,ω,w)
		k2 = f_dir_cycle_k(θ+k1*h/2,ω,w)
		k3 = f_dir_cycle_k(θ+k2*h/2,ω,w)
		k4 = f_dir_cycle_k(θ+k3*h,ω,w)

		dθ = (k1+2*k2+2*k3+k4)./6
		θ += h*dθ
		corr = norm(dθ,Inf)
	end

	return θ
end

function dir_cycle_kuramoto(θ0::Vector{Float64}, ω::Vector{Float64}, h::Float64=.01, max_iter::Int64=1000, thr::Float64=1e-6)
	return dir_cycle_kuramoto(θ0,ω,ones(length(θ0)),h,max_iter,thr)
end

function f_dir_cycle_k(θ::Vector{Float64}, ω::Vector{Float64}, w::Vector{Float64})
	n = length(θ)

	return ω - w.*sin.(θ - [θ[2:end];θ[1]])
end

function f_dir_cycle_k(θ::Vector{Float64}, ω::Vector{Float64})
	return f_dir_cycle_k(θ,ω,ones(length(θ)))
end


