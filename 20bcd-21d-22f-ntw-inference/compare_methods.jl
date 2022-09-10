using PyPlot, Statistics

include("rand_graphs.jl")
include("line_dist.jl")

n = 20

n_iter = 10
n_range = 10
h = .01
τmin = 1.1*h
τmax = 3.
τs = LinRange(τmin,τmax,n_range)

M1 = Matrix{Float64}(undef,n_iter,0)
M2 = Matrix{Float64}(undef,n_iter,0)

for τ in τs
	m1 = Vector{Float64}()
	m2 = Vector{Float64}()
	for i in 1:n_iter
		global n,h

		B = Matrix(my_er(n,.2))
		n,m = size(B)
		L = B*B'
		ω = .3*rand(n)
		ω .-= mean(ω)

		a0 = -1.

	#	L *= 10
	#	a0 *= 10

		Ti = 10*h
		Tf = Ti + τ

		x0 = pinv(L)*ω
		ψ0 = L*x0

	 	e = rand(1:m)
		ij = (findmax(B[:,e])[2],findmin(B[:,e])[2])
	 #=	
		xs = linear_step(L,ω,x0,ij,a0,Ti,Tf,true,1000,-1.,h) 
	# =#
	# #=
		x0 = vec(kuramoto_step(L,ω,x0,e,0.,0.,0.,false,10000,1e-5,h))
		xs = kuramoto_step(L,ω,x0,e,a0,Ti,Tf,true,1000,-1.,h)
	# =#
		ψs = L*xs
		ψ0 = L*x0

		dx = Matrix{Float64}(undef,n,0)
		dψ = Matrix{Float64}(undef,n,0)
		for t in 1:size(xs)[2]
			dx = [dx xs[:,t]-x0]
			dψ = [dψ ψs[:,t]-ψ0]
		end

		 #= # 1-norm
		X = vec(sum(abs.(dx),dims=2))
		Ψ = vec(sum(abs.(dψ),dims=2))
		# =#
		# #= ∞-norm
		X = [maximum(abs.(dx[i,:])) for i in 1:n]
		Ψ = [maximum(abs.(dψ[i,:])) for i in 1:n]
		# =#

		ix = findmax(X)[2]
		jx = findmax([X[1:ix-1];0.;X[ix+1:n]])[2]
		iψ = findmax(Ψ)[2]
		jψ = findmax([Ψ[1:iψ-1];0.;Ψ[iψ+1:n]])[2]

		e0 = sort([ij[1],ij[2]])
		ex = sort([ix,jx])
		eψ = sort([iψ,jψ])

		ux = union(e0,ex)
		if length(ux) == 4
			push!(m1,0.)
		elseif length(ux) == 3
			push!(m1,.5)
		elseif length(ux) == 2
			push!(m1,1.)
		end

		uψ = union(e0,eψ)
		if length(uψ) == 4
			push!(m2,0.)
		elseif length(uψ) == 3
			push!(m2,.5)
		elseif length(uψ) == 2
			push!(m2,1.)
		end
	end

	global M1 = [M1 m1]
	global M2 = [M2 m2]
end

PyPlot.plot(τs,vec(sum(M1,dims=1))./n_iter,"--",color="C0")
PyPlot.plot(τs,vec(sum((M1 .== 1.)./n_iter,dims=1)),"-",color="C0")
PyPlot.plot(τs,vec(sum(M2,dims=1))./n_iter,"--",color="C1")
PyPlot.plot(τs,vec(sum((M2 .== 1.)./n_iter,dims=1)),"-",color="C1")

xlabel("τ")
ylabel("ρ")
title("n = $n, n_iter = $n_iter")





