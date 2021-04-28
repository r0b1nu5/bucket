using DelimitedFiles, LinearAlgebra, LightGraphs, OrdinaryDiffEq, NetworkDynamics

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

	dθs = Array{Float64,2}(undef,n,0)

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
			
			writedlm("temp_data/temp_θ_$c.csv",θs[:,1:end-1],',')
			θ = θs[:,end]

			writedlm("temp_data/temp_dθ_$c.csv",dθs[:,1:end],',')
			dθs = Array{Float64,2}(undef,n,0)
		end

		k1 = ω - B12*WW*(sin.(BB'*θ .- α) .+ sin(α))
		k2 = ω - B12*WW*(sin.(BB'*(θ + h/2*k1) .- α) .+ sin(α))
		k3 = ω - B12*WW*(sin.(BB'*(θ + h/2*k2) .- α) .+ sin(α))
		k4 = ω - B12*WW*(sin.(BB'*(θ + h*k3) .- α) .+ sin(α))

		dθ = (k1 + 2*k2 + 2*k3 + k4)/6

		θ += h*dθ

		θs = [θs θ]
		dθs = [dθs dθ]

		err = maximum(dθ)-minimum(dθ)
	end

	Θs = Array{Float64,2}(undef,n,0)
	dΘs = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Θs = [Θs readdlm("temp_data/temp_θ_$i.csv",',')]
		rm("temp_data/temp_θ_$i.csv")
		dΘs = [dΘs readdlm("temp_data/temp_dθ_$i.csv",',')]
		rm("temp_data/temp_dθ_$i.csv")
	end
	Θs = [Θs θs]
	dΘs = [dΘs dθs]

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

	dθs = Array{Float64,2}(undef,n,0)

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
			
			writedlm("temp_data/temp_θ_$c.csv",θs[:,1:end-1],',')
			θs = θs[:,end]

			writedlm("temp_data/temp_dθ_$c.csv",dθs[:,1:end],',')
			dθs = Array{Float64,2}(undef,n,0)
		end

		k1 = ω - B12*WW*(sin.(BB'*θ .- α) .+ sin(α))
		k2 = ω - B12*WW*(sin.(BB'*(θ + h/2*k1) .- α) .+ sin(α))
		k3 = ω - B12*WW*(sin.(BB'*(θ + h/2*k2) .- α) .+ sin(α))
		k4 = ω - B12*WW*(sin.(BB'*(θ + h*k3) .- α) .+ sin(α))

		dθ = (k1 + 2*k2 + 2*k3 + k4)/6

		θ += h*dθ

		θs = [θs θ]
		dθs = [dθs dθ]

		err = maximum(dθ)-minimum(dθ)
	end

	Θs = Array{Float64,2}(undef,n,0)
	dΘs = Array{Float64,2}(undef,n,0)
	for i in 1:c
		Θs = [Θs readdlm("temp_data/temp_θ_$i.csv",',')]
		rm("temp_data/temp_θ_$i.csv")
		dΘs = [dΘs readdlm("temp_data/temp_dθ_$i.csv",',')]
		rm("temp_data/temp_dθ_$i.csv")
	end
	Θs = [Θs θs]
	dΘs = [dΘs dθs]

	return Θs,dΘs,err,iter
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

function ksakaguchi_ND(L::Array{Float64,2}, ω::Array{Float64,1}, θ0::Array{Float64,1}, α::Float64, t_span::Tuple{Float64,Float64}, dts::Tuple{Float64,Float64,Float64}, algo::Any=RK4(), verb::Bool=false)
	n = length(ω)

	A = -L + diagm(0 => diag(L))

	G = SimpleDiGraph(A)

	ks_sol = solve(load_ksakaguchi(G,ω,θ0,t_span,dts,α,1.),algo)

	ts = ks_sol.t
	T = length(ts)

	θ = Array{Float64,2}(undef,n,0)
	for i in 1:T
		θ = [θ ks_sol.u[i]]
	end

	return ts, θ
end


function ksakaguchi_ND(L::SparseMatrixCSC{Float64,Int}, ω::Array{Float64,1}, θ0::Array{Float64,1}, α::Float64, t_span::Tuple{Float64,Float64}, dts::Tuple{Float64,Float64,Float64}, algo::Any=RK4(), verb::Bool=false)
	n = length(ω)

	A = -L + spdiagm(0 => diag(L))

	G = SimpleDiGraph(A)

	ks_sol = solve(load_ksakaguchi(G,ω,θ0,t_span,dts,α,1.),algo)

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

function load_ksakaguchi(G::SimpleDiGraph{Int64}, ω::Array{Float64,1}, θ0::Array{Float64,1}, t_span::Tuple{Float64,Float64}, dts::Tuple{Float64,Float64,Float64}, α::Float64=.1, K::Float64=1.)
	vertex! = ODEVertex(f! = ks_vertex!, dim=1)
	edge! = StaticEdge(f! = ks_edge!, dim=1, coupling = :directed)

	ksakaguchi! = network_dynamics(vertex!, edge!, G)

	vpar = ω
	epar = (K,α)
	p = (vpar,epar)

	return ODEProblem(ksakaguchi!, θ0, t_span, p, dt=dts[1], dtmax=dts[2], dtmin=dts[3])
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


function freq_width(L::Array{Float64,2}, ω0::Array{Float64,1}, θ0::Array{Float64,1}, α::Float64, C::Array{Int64,1}, verb::Bool=false, res::Float64=.0005)
	if norm(ω0) < 1e-8
		@info "The frequency vector is close to zero."
	end

	n = length(θ0)
	
	ω = ω0 .- mean(ω0)
	ω /= norm(ω)

	x = ksakaguchi(L,zeros(n),θ0,α,true,false,.01,1e-6)
#	ts,x = ksakaguchi_ND(L,zeros(n),θ0,α,(0.,10.))
	θ1 = x[1][:,end]
	θ = copy(θ1)
	q0 = winding(θ,C)

	β = 0.
	dβ = 1.

	while dβ > res
		if verb
			@info "β = $β, dβ = $dβ"
		end

		q = q0
		it = 0
		while q == q0 && it < 100000
			β += dβ
			x = ksakaguchi(L,β*ω,θ,α,true,false,.01,1e-6)
#			ts,x = ksakaguchi_ND(L,β*ω,θ,α,(0.,10.))
			q = winding(x[1][:,end],Array(1:n))
			it = x[4]
			
			if verb
				@info "q = $q, it = $it"
			end
		end
		θ = x[1][:,1]
		β -= dβ
		dβ /= 10
	end

	βmax = copy(β)
	fmax = mean(x[2][:,end])

	β = 0.
	dβ = 1.

	while dβ > res
		if verb
			@info "β = $β, dβ = $dβ"
		end

		q = q0
		it = 0
		while  q == q0 && it < 100000
			β -= dβ
			x = ksakaguchi(L,β*ω,θ,α,true,false,.01,1e-6)
#			ts,x = ksakaguchi_ND(L,β*ω,θ,α,(0.,10.))
			q = winding(x[1][:,end],Array(1:n))
			it = x[4]

			if verb
				@info "q = $q, it = $it"
			end
		end
		θ = x[1][:,1]
		β += dβ
		dβ /= 10
	end

	βmin = copy(β)
	fmin = mean(x[2][:,end])

	return βmin,βmax,fmin,max
end

# Loads the coupling function for the Kuramoto-Sakaguchi model.

function h(x::Float64, α::Float64=.1)
	return sin(x - α) + sin(α)
end

function h(x::Union{Array{Float64,1},LinRange{Float64}}, α::Float64=.1)
	return [h(x[i]) for i in 1:length(x)]
end

function hi(f::Float64, α::Float64=.1)
	return asin(f - sin(α)) + α
end

function hi(f::Union{Array{Float64,1},LinRange{Float64}}, α::Float64=.1)
	return [hi(f[i]) for i in 1:length(f)]
end

function H(f::Float64, α::Float64=.1)
	return h(-hi(f,α),α)
end

function H(f::Union{Array{Float64,1},LinRange{Float64}}, α::Float64=.1)
	return [H(f[i]) for i in 1:length(f)]
end


