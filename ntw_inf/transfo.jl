using PyPlot, LinearAlgebra, DelimitedFiles

include("../L2B.jl")

# Runs Kuramoto on L with coupling on line ij subject to oscillations

function kuramoto_transfo(L::Array{Float64,2}, P::Array{Float64,1}, th0::Array{Float64,1}, ij::Tuple{Int64,Int64}, osc::Tuple{Float64,Float64,Float64}, T::Int64, history::Bool=true, eps::Float64=1e-6, h::Float64=.1)
	n = size(L)[1]
	
	i,j = ij
	B,w,kij = L2B_ij(L,i,j)
	Bt = Array(transpose(B))
	W = diagm(0 => w)

	a0,w0,p0 = osc
	eij = zeros(n)
	eij[[i,j]] = [1,-1]

	ths = Array{Float64,2}(undef,n,0)
	dths = Array{Float64,2}(undef,n,0)

	error = 1000.
	iter = 0
	
	th1 = copy(th0)
	th2 = copy(th1)

	Wt = copy(W)

	it = 0
	while error > eps
		it += 1
		th1 = copy(th2)

		k1 = P - B*W*sin.(Bt*th1)
		k2 = P - B*W*sin.(Bt*(th1+h/2*k1))
		k3 = P - B*W*sin.(Bt*(th1+h/2*k2))
		k4 = P - B*W*sin.(Bt*(th1+h*k3))

		dth = (k1+2*k2+2*k3+k4)/6
		th2 = th1 + h*dth

		error = maximum(abs.(dth))

		if it%1000 == 0
			@info "$error"
		end
	end

	while iter < T
		iter += 1
		if iter%=1000 == 0
			@info "$iter"
		end

		h1 = copy(th2)

		xi = a0*sin(w0*h*iter + p0)
		Wt[kij,kij] = W[kij,kij]*(1 + xi)
		
		k1 = P - B*Wt*sin.(Bt*th1)
		k2 = P - B*Wt*sin.(Bt*(th1+h/2*k1))
		k3 = P - B*Wt*sin.(Bt*(th1+h/2*k2))
		k4 = P - B*Wt*sin.(Bt*(th1+h*k3))

		dth = (k1+2*k2+2*k3+k4)/6
		th2 = th1 + h*dth

		ths = [ths th2]
		dths = [dths dth]
	end

	return ths,dths
end


# Computes the transformed time series on the observed nodes. 
# Nodes in Ic are Kron-reduced.

function psi(L::Array{Float64,2}, Ic::Array{Int64,1}, ths::Array{Float64,2})
	n = size(L)[1]

	Ig = setdiff(Array(1:n),Ic)

	Lgg = L[Ig,Ig]
	Lgc = L[Ig,Ic]
	Lcg = L[Ic,Ig]
	Lcc = L[Ic,Ic]

	Lr = Lgg - Lgc*inv(Symmetric(Lcc))*Lcg

	ps = Lr*ths[Ig,:]

	return ps
end











