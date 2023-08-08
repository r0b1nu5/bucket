using PyPlot

include("tools.jl")

function dir_kuramoto(L::Matrix{Float64}, θ0::Vector{Float64}, ω::Vector{Float64}, h = .01, max_iter::Int64=1000, thr::Float64=1e-6)
	Bout,Bin,w = L2B(L)

	θ = θ0

	t = 0
	corr = 1000.
	while t < max_iter && (corr > thr || t < 10)
		t += 1
#		if t%100 == 0
#			@info "$t/$max_iter"
#		end

		k1 = f_dk(θ,ω,Bout,Bin,w)
		k2 = f_dk(θ+k1*h/2,ω,Bout,Bin,w)
		k3 = f_dk(θ+k2*h/2,ω,Bout,Bin,w)
		k4 = f_dk(θ+k3*h,ω,Bout,Bin,w)

		dθ = (k1+2*k2+2*k3+k4)./6
		θ += h*dθ
		corr = norm(dθ,Inf)
	end

	return θ
end

function f_dk(θ::Vector{Float64}, ω::Vector{Float64}, Bout::Matrix{Float64}, Bin::Matrix{Float64}, w::Vector{Float64})
	return ω - Bout*diagm(0 => w)*sin.((Bout-Bin)'*θ)
end

function f_dk(θ::Vector{Float64}, ω::Vector{Float64}, Bout::Matrix{Float64}, Bin::Matrix{Float64})
	return f_dk(θ,ω,Bout,Bin,ones(size(Bout)[2]))
end


	



