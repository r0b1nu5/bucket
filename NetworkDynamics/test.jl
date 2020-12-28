Ns = [5,10,20,50,100,200,500,1000,2000,5000,10000]

T1 = Array{Float64,1}()
T2 = Array{Float64,1}()

for n in Ns
	local g = watts_strogatz(n, 3, 0.2)
	local A = adjacency_matrix(g)
	local G = SimpleDiGraph(A)
	local L = Float64.(laplacian_matrix(G))
	
	local om = .2*rand(n) .- .1
	local th0 = 2pi*rand(n)
	local Ti = 0.
	local Tf = 40.
	local tspan = (Ti,Tf)

	t1 = time()
	local ks_sol = solve(load_ksakaguchi(G,om,th0,tspan,.5,1.),Tsit5())
	t2 = time()
@info "$(now()) -- ND.jl, n=$n done."

	t3 = time()
	local xxx = ksakaguchi(L,om,th0,a,true,true,.01,-1.,4000)
	t4 = time()
@info "$(now()) -- Own, n=$n done."
	
	push!(T1,t2-t1)
	push!(T2,t4-t3)
end

PyPlot.plot(Ns,T1,Ns,T2)





