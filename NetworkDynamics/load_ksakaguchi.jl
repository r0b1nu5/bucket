using LightGraphs, OrdinaryDiffEq, NetworkDynamics

#n = 200
#g = watts_strogatz(n, 3, 0.2)
#A = adjacency_matrix(g)
#G = SimpleDiGraph(A)

function load_ksakaguchi(G::SimpleDiGraph{Int64}, ω::Array{Float64,1}, θ0::Array{Float64,1}, t_span::Tuple{Float64,Float64}, α::Float64=.1, K::Float64=1.)
	vertex! = ODEVertex(f! = ks_vertex!, dim=1)
	edge! = StaticEdge(f! = ks_edge!, dim=1, coupling = :directed)

	ksakaguchi! = network_dynamics(vertex!, edge!, G)

	vpar = ω
	epar = (K,α)
	p = (vpar,epar)

	return ODEProblem(ksakaguchi!, θ0, t_span, p)
end

function ks_edge!(e, θ_s, θ_d, epar, t)
	e .= -epar[1].*sin.(-θ_s .+ θ_d .- epar[2])
end

function ks_vertex!(dθ, θ, e_s, vpar, t)
	dth .= vpar
	for edge in e_s
		dθ .+= edge
	end
end






