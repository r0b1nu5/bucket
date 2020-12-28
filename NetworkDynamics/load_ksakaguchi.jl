using LightGraphs, OrdinaryDiffEq, NetworkDynamics

#n = 200
#g = watts_strogatz(n, 3, 0.2)
#A = adjacency_matrix(g)
#G = SimpleDiGraph(A)

function load_ksakaguchi(G::SimpleDiGraph{Int64}, om::Array{Float64,1}, th0::Array{Float64,1}, tspan::Tuple{Float64,Float64}, a::Float64=.1, K::Float64=1.)
	vertex! = ODEVertex(f! = ks_vertex!, dim=1)
	edge! = StaticEdge(f! = ks_edge!, dim=1)

	ksakaguchi! = network_dynamics(vertex!, edge!, G)

	vpar = om
	epar = (K,a)
	p = (vpar,epar)

	return ODEProblem(ksakaguchi!, th0, tspan, p)
end

function ks_edge!(e, th_s, th_d, epar, t)
	e .= -epar[1].*sin.(th_s .- th_d .- epar[2])
end

function ks_vertex!(dth, th, e_s, e_d, vpar, t)
	dth .= vpar
	for edge in e_s
		dth .+= edge
	end
end






