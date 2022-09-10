using PyPlot, LinearAlgebra, DelimitedFiles

include("tools.jl")

#TODO Using sine coupling...


# !!! Directed graph !!!
function coev(B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, x0::Vector{Float64}, a0::Vector{Float64}, save::Bool=true, δmin::Float64=0., δmax::Float64=1., γ::Float64=1., σ::Float64=100., tol::Float64=1e-6, maxiter::Int64=10000, h::Float64=.001)
	n,m = size(B)
	Bout = B.*(B .> 1e-2)

	Δ = δmax - δmin

	xs = copy(x0)
	as = copy(a0)

	iter = 0
	c = 0
	err = 1000.

	while iter < maxiter && err > tol
		iter += 1

		if iter%1000 == 0
			@info "iter: $iter, err: $err"
			c += 1
			writedlm("temp/xs_$c.csv",xs[:,1:end-1],',')
			xs = xs[:,end]
			#writedlm("temp/as_$c.csv",as[:,1:end-1],',')
			as = as[:,end]
		end

		x = xs[:,end]
		dx = B'*x
		a = as[:,end]

		kx1 = kx(x,a,B,Bout)
		ka1 = ka(x,a,B,σ,δmin,δmax)/γ
		kx2 = kx(x+h/2*kx1,a+h/2*ka1,B,Bout)
		ka2 = ka(x+h/2*kx1,a+h/2*ka1,B,σ,δmin,δmax)/γ
		kx3 = kx(x+h/2*kx2,a+h/2*ka2,B,Bout)
		ka3 = ka(x+h/2*kx2,a+h/2*ka2,B,σ,δmin,δmax)/γ
		kx4 = kx(x+h*kx3,a+h*ka3,B,Bout)
		ka4 = ka(x+h*kx3,a+h*ka3,B,σ,δmin,δmax)/γ

		dx = (kx1 + 2*kx2 + 2*kx3 + kx4)/6
		da = (ka1 + 2*ka2 + 2*ka3 + ka4)/6

		err = max(maximum(abs.(dx)),maximum(abs.(da)))

		xs = [xs x+h*dx]
		as = [as a+h*da]
	end

	Xs = Matrix{Float64}(undef,n,0)
	As = Matrix{Float64}(undef,m,0)

	for i in 1:c
		Xs = [Xs readdlm("temp/xs_$i.csv",',')]
		rm("temp/xs_$i.csv")
		#As = [As readdlm("temp/as_$i.csv",',')]
		#rm("temp/as_$i.csv")
	end
	Xs = [Xs xs]
	As = [As as]

	return Xs, As, iter, err
end

function kx(x::Vector{Float64}, a::Vector{Float64}, B::Matrix{Float64}, Bout::Matrix{Float64})
	return -Bout*diagm(0 => a)*B'*x
end

function kx(x::Vector{Float64}, a::Vector{Float64}, B::SparseMatrixCSC{Float64,Int64}, Bout::SparseMatrixCSC{Float64,Int64})
	return -Bout*spdiagm(0 => a)*(x'*B)'
end

function ka(x::Vector{Float64}, a::Vector{Float64}, B::Union{Matrix{Float64},SparseMatrixCSC{Float64,Int64}}, σ::Float64=1., δmin::Float64=0., δmax::Float64=1.)
	return vec((δmax-δmin)*exp.(-σ*(x'*B).^2)) .+ δmin - a
end





