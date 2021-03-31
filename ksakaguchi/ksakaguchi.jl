using DelimitedFiles, LinearAlgebra, LightGraphs, OrdinaryDiffEq, NetworkDynamics

include("L2B.jl")

function ksakaguchi(L::Array{Float64,2}, om::Array{Float64,1}, th0::Array{Float64,1}, a::Float64, save_history::Bool=false, verb::Bool=false, h::Float64=.01, thres::Float64=1e-5, max_iter::Int64=100000)
	B,w = L2B(L)
	W = diagm(0 => w)
	n,m = size(B)
	B1 = B.*(B .> 0)
	B2 = -B.*(B .< 0)
	B12 = [B1 B2]
	BB = [B -B]
	WW = [W zeros(m,m);zeros(m,m) W]
	
	th = th0
	ths = th0

	dhs = Array{Float64,2}(undef,n,0)

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
			
			writedlm("temp_data/temp_th_$c.csv",ths[:,1:end-1],',')
			ths = ths[:,end]

			writedlm("temp_data/temp_dh_$c.csv",dhs[:,1:end],',')
			dhs = Array{Float64,2}(undef,n,0)
		end

		k1 = om - B12*WW*sin.(BB'*th .- a)
		k2 = om - B12*WW*sin.(BB'*(th + h/2*k1) .- a)
		k3 = om - B12*WW*sin.(BB'*(th + h/2*k2) .- a)
		k4 = om - B12*WW*sin.(BB'*(th + h*k3) .- a)

		dh = (k1 + 2*k2 + 2*k3 + k4)/6

		th += h*dh

		ths = [ths th]
		dhs = [dhs dh]

		err = maximum(dh)-minimum(dh)
	end

	Ths = Array{Float64,2}(undef,n,0)
	Dhs = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("temp_data/temp_th_$i.csv",',')]
		rm("temp_data/temp_th_$i.csv")
		Dhs = [Dhs readdlm("temp_data/temp_dh_$i.csv",',')]
		rm("temp_data/temp_dh_$i.csv")
	end
	Ths = [Ths ths]
	Dhs = [Dhs dhs]

	return Ths,Dhs,err,iter
end

function ksakaguchi(L::SparseMatrixCSC{Float64,Int}, om::Array{Float64,1}, th0::Array{Float64,1}, a::Float64, save_history::Bool=false, verb::Bool=false, h::Float64=.01, thres::Float64=1e-5, max_iter::Int64=100000)
	B,w,Bt = L2B(L)
	W = spdiagm(0 => w)
	n,m = size(B)
	B1 = B.*(B .> 0)
	B2 = -B.*(B .< 0)
	B12 = [B1 B2]
	BB = [B -B]
	WW = [W zeros(m,m);zeros(m,m) W]
	
	th = th0
	ths = th0

	dhs = Array{Float64,2}(undef,n,0)

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
			
			writedlm("temp_data/temp_th_$c.csv",ths[:,1:end-1],',')
			ths = ths[:,end]

			writedlm("temp_data/temp_dh_$c.csv",dhs[:,1:end],',')
			dhs = Array{Float64,2}(undef,n,0)
		end

		k1 = om - B12*WW*sin.(BB'*th .- a)
		k2 = om - B12*WW*sin.(BB'*(th + h/2*k1) .- a)
		k3 = om - B12*WW*sin.(BB'*(th + h/2*k2) .- a)
		k4 = om - B12*WW*sin.(BB'*(th + h*k3) .- a)

		dh = (k1 + 2*k2 + 2*k3 + k4)/6

		th += h*dh

		ths = [ths th]
		dhs = [dhs dh]

		err = maximum(dh)-minimum(dh)
	end

	Ths = Array{Float64,2}(undef,n,0)
	Dhs = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Ths = [Ths readdlm("temp_data/temp_th_$i.csv",',')]
		rm("temp_data/temp_th_$i.csv")
		Dhs = [Dhs readdlm("temp_data/temp_dh_$i.csv",',')]
		rm("temp_data/temp_dh_$i.csv")
	end
	Ths = [Ths ths]
	Dhs = [Dhs dhs]

	return Ths,Dhs,err,iter
end


function ksakaguchi_ND(L::Array{Float64,2}, ω::Array{Float64,1}, θ0::Array{Float64,1}, α::Float64, t_span::Tuple{Float64,Float64}, algo::Any=RK4(), verb::Bool=false)
	n = length(ω)

	A = -L + diagm(0 => diag(L))

	G = SimpleDiGraph(A)

	ks_sol = solve(load_ksakaguchi(G,ω,θ0,t_span,α,1.),algo)

	ts = ks_sol.t
	T = length(ts)

	θ = Array{Float64,2}(undef,n,0)
	for i in 1:T
		θ = [θ ks_sol.u[i]]
	end

	return ts, θ
end


function ksakaguchi_ND(L::SparseMatrixCSC{Float64,Int}, ω::Array{Float64,1}, θ0::Array{Float64,1}, α::Float64, t_span::Tuple{Float64,Float64}, algo::Any=RK4(), verb::Bool=false)
	n = length(ω)

	A = -L + spdiagm(0 => diag(L))

	G = SimpleDiGraph(A)

	ks_sol = solve(load_ksakaguchi(G,ω,θ0,t_span,α,1.),algo)

	ts = ks_sol.t
	T = length(ts)

	θ = Array{Float64,2}(undef,n,0)
	for i in 1:T
		θ = [θ ks_sol.u[i]]
	end

	return ts, θ
end

function load_ksakaguchi(G::SimpleDiGraph{Int64}, ω::Array{Float64,1}, θ0::Array{Float64,1}, t_span::Tuple{Float64,Float64}, α::Float64=.1, K::Float64=1.)
	vertex! = ODEVertex(f! = ks_vertex!, dim=1)
	edge! = StaticEdge(f! = ks_edge!, dim=1, coupling = :directed)

	ksakaguchi! = network_dynamics(vertex!, edge!, G)

	vpar = ω
	epar = (K,α)
	p = (vpar,epar)

	return ODEProblem(ksakaguchi!, θ0, t_span, p)
end

function ks_edge!(e, θ_s, θ_t, epar, t)
	e .= -epar[1].*sin.(θ_t .- θ_s .- epar[2])
end

function ks_vertex!(dθ, θ, e_s, vpar, t)
	dθ .= vpar
	for edge in e_s
		dθ .+= edge
	end
end








