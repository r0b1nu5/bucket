using DelimitedFiles, LinearAlgebra, PyPlot
using LightGraphs, OrdinaryDiffEq, NetworkDynamics

include("tools.jl")


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




