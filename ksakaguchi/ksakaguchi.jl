using DelimitedFiles, LinearAlgebra, PyPlot

include("tools.jl")

function ksakaguchi(L::Array{Float64,2}, ω::Array{Float64,1}, θ0::Array{Float64,1}, α::Float64, save_history::Bool=false, verb::Bool=false, h::Float64=.01, thres::Float64=1e-5, max_iter::Int64=100000)
	B,w = L2B(L)
	W = diagm(0 => w)
	n,m = size(B)
	B1 = B.*(B .> 0)
	B2 = -B.*(B .< 0)
	B12 = [B1 B2]
	BB = [B -B]
	WW = [W zeros(m,m);zeros(m,m) W]
	
	θ = θ0
	θs = θ0

	dθ = zeros(n)
	dθs = Array{Float64,2}(undef,n,0)

	err = 1000.
	iter = 0
	c = 0

	while err > thres && iter < max_iter
		iter += 1

		if iter%1000 == 0
			if verb
				@info "iter: $iter, err = $(round(err,digits=5))"
			end

			if save_history
				c += 1
				
				writedlm("temp_data/temp_θ_$c.csv",θs[:,1:end-1],',')
				θs = θs[:,end]
	
				writedlm("temp_data/temp_dθ_$c.csv",dθs[:,1:end],',')
				dθs = Array{Float64,2}(undef,n,0)
			end
		end

		k1 = ω - B12*WW*(sin.(BB'*θ .- α) .+ sin(α))
		k2 = ω - B12*WW*(sin.(BB'*(θ + h/2*k1) .- α) .+ sin(α))
		k3 = ω - B12*WW*(sin.(BB'*(θ + h/2*k2) .- α) .+ sin(α))
		k4 = ω - B12*WW*(sin.(BB'*(θ + h*k3) .- α) .+ sin(α))

		dθ = (k1 + 2*k2 + 2*k3 + k4)/6

		θ += h*dθ

		if save_history
			θs = [θs θ]
			dθs = [dθs dθ]
		end

		err = maximum(dθ)-minimum(dθ)
	end

	Θs = Array{Float64,2}(undef,n,0)
	dΘs = Array{Float64,2}(undef,n,0)
	if save_history
		for i in 1:c
			Θs = [Θs readdlm("temp_data/temp_θ_$i.csv",',')]
			rm("temp_data/temp_θ_$i.csv")
			dΘs = [dΘs readdlm("temp_data/temp_dθ_$i.csv",',')]
			rm("temp_data/temp_dθ_$i.csv")
		end
		Θs = [Θs θs]
		dΘs = [dΘs dθs]
	else
		Θs = [Θs θ]
		dΘs = [dΘs dθ]
	end

	return Θs,dΘs,err,iter
end

function ksakaguchi(L::SparseMatrixCSC{Float64,Int}, ω::Array{Float64,1}, θ0::Array{Float64,1}, α::Float64, save_history::Bool=false, verb::Bool=false, h::Float64=.01, thres::Float64=1e-5, max_iter::Int64=100000)
	B,w,Bt = L2B(L)
	W = spdiagm(0 => w)
	n,m = size(B)
	B1 = B.*(B .> 0)
	B2 = -B.*(B .< 0)
	B12 = [B1 B2]
	BB = [B -B]
	WW = [W zeros(m,m);zeros(m,m) W]
	
	θ = θ0
	θs = θ0

	dθ = zeros(n)
	dθs = Array{Float64,2}(undef,n,0)

	err = 1000.
	iter = 0
	c = 0

	while err > thres && iter < max_iter
		iter += 1

		if iter%1000 == 0
			if verb
				@info "iter: $iter, err = $(round(err,digits=5))"
			end
			
			if save_history
				c += 1
				writedlm("temp_data/temp_θ_$c.csv",θs[:,1:end-1],',')
				θs = θs[:,end]
	
				writedlm("temp_data/temp_dθ_$c.csv",dθs[:,1:end],',')
				dθs = Array{Float64,2}(undef,n,0)
			end
		end

		k1 = ω - B12*WW*(sin.(BB'*θ .- α) .+ sin(α))
		k2 = ω - B12*WW*(sin.(BB'*(θ + h/2*k1) .- α) .+ sin(α))
		k3 = ω - B12*WW*(sin.(BB'*(θ + h/2*k2) .- α) .+ sin(α))
		k4 = ω - B12*WW*(sin.(BB'*(θ + h*k3) .- α) .+ sin(α))

		dθ = (k1 + 2*k2 + 2*k3 + k4)/6

		θ += h*dθ

		if save_history
			θs = [θs θ]
			dθs = [dθs dθ]
		end

		err = maximum(dθ)-minimum(dθ)
	end

	Θs = Array{Float64,2}(undef,n,0)
	dΘs = Array{Float64,2}(undef,n,0)
	if save_history
		for i in 1:c
			Θs = [Θs readdlm("temp_data/temp_θ_$i.csv",',')]
			rm("temp_data/temp_θ_$i.csv")
			dΘs = [dΘs readdlm("temp_data/temp_dθ_$i.csv",',')]
			rm("temp_data/temp_dθ_$i.csv")
		end
		Θs = [Θs θs]
		dΘs = [dΘs dθs]
	else
		Θs = [Θs θ]
		dΘs = [dΘs dθ]
	end

	return Θs,dΘs,err,iter
end


function ks_flows(θ::Array{Float64,1}, L::Array{Float64,2}, α::Float64)
	B,w = L2B(L)
	W = diagm(0 => w)
	n,m = size(B)
	B1 = B.*(B .> 0)
	B2 = -B.*(B .< 0)
	B12 = [B1 B2]
	BB = [B -B]
	WW = [W zeros(m,m);zeros(m,m) W]
	
	return WW*(sin.(BB'*θ .- α) .+ sin(α))
end


# Loads the coupling function for the Kuramoto-Sakaguchi model.

#=
function h(x::Float64, α::Float64=.1)
	return sin(x - α) + sin(α)
end

function h(x::Union{Array{Float64,1},LinRange{Float64}}, α::Float64=.1)
	return [h(x[i],α) for i in 1:length(x)]
end

function h(x::Array{Float64,2}, α::Float64=.1)
	hh = Array{Float64,2}(undef,size(x)[1],0)
	for j in 1:size(x)[2]
		hh = [hh h(x[:,j],α)]
	end

	return hh
end

function hi(f::Float64, α::Float64=.1)
	return asin(f - sin(α)) + α
end

function hi(f::Union{Array{Float64,1},LinRange{Float64}}, α::Float64=.1)
	return [hi(f[i],α) for i in 1:length(f)]
end

function hi(f::Array{Float64,2}, α::Float64=.1)
	hh = Array{Float64,2}(undef,size(f)[1],0)
	for j in 1:size(f)[2]
		hh = [hh hi(f[:,j],α)]
	end

	return hh
end

function H(f::Float64, α::Float64=.1)
	return h(-hi(f,α),α)
end

function H(f::Union{Array{Float64,1},LinRange{Float64}}, α::Float64=.1)
	return [H(f[i]) for i in 1:length(f)]
end
=#

function ksakaguchi_2nd(L::Array{Float64,2}, ω::Array{Float64,1}, d::Array{Float64,1}, θ0::Array{Float64,1}, θd0::Array{Float64,1}, α::Float64, save_history::Bool=false, verb::Bool=false, h::Float64=.01, thres::Float64=1e-5, max_iter::Int64=100000)
	B,w = L2B(L)
	W = diagm(0 => w)
	n,m = size(B)
	B1 = B.*(B .> 0)
	B2 = -B.*(B .< 0)
	B12 = [B1 B2]
	BB = [B -B]
	WW = [W zeros(m,m);zeros(m,m) W]
	D = diagm(0 => d)
	
	x = [θ0;θd0]
	xs = x

	dxs = Array{Float64,2}(undef,2*n,0)

	err = 1000.
	iter = 0
	c = 0

	while err > thres && iter < max_iter
		iter += 1

		if iter%1000 == 0
			c += 1
			if verb
				@info "iter: $iter, err = $(round(err,digits=5))"
			end
	
			writedlm("temp_data/temp_x_$c.csv",xs[:,1:end-1],',')
			xs = xs[:,end]

			writedlm("temp_data/temp_dx_$c.csv",dxs[:,1:end],',')
			dxs = Array{Float64,2}(undef,2*n,0)
		end

		θ = x[1:n]
		θd = x[(n+1):(2*n)]

		k11 = θd
		k21 = ω - D*θd - B12*WW*(sin.(BB'*θ .- α) .+ sin(α))
		k12 = θd + h/2*k21
		k22 = ω - D*(θd + h/2*k21) - B12*WW*(sin.(BB'*(θ + h/2*k11) .- α) .+ sin(α))
		k13 = θd + h/2*k22
		k23 = ω - D*(θd + h/2*k22) - B12*WW*(sin.(BB'*(θ + h/2*k12) .- α) .+ sin(α))
		k14 = θd + h*k23
		k24 = ω - D*(θd + h*k23) - B12*WW*(sin.(BB'*(θ + h*k13) .- α) .+ sin(α))

		dθ = (k11 + 2*k12 + 2*k13 + k14)/6
		dθd = (k21 + 2*k22 + 2*k23 + k24)/6
		dx = [dθ;dθd]
		x += h*dx

		xs = [xs x]
		dxs = [dxs dx]

		err = maximum(dx)-minimum(dx)
	end

	Xs = Array{Float64,2}(undef,2*n,0)
	dXs = Array{Float64,2}(undef,2*n,0)
	for i in 1:c
		Xs = [Xs readdlm("temp_data/temp_x_$i.csv",',')]
		rm("temp_data/temp_x_$i.csv")
		dXs = [dXs readdlm("temp_data/temp_dx_$i.csv",',')]
		rm("temp_data/temp_dx_$i.csv")
	end
	Xs = [Xs xs]
	dXs = [dXs dxs]

	return Xs,dXs,err,iter
end


