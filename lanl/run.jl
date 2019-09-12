using Distributed

n_thr = 45

if nworkers() < n_thr
	addprocs(n_thr - nworkers())
end

@everywhere include("final.jl")

@everywhere function locate(id::Int64, n::Int64, a0::Float64, f0::Float64, p0::Float64, T::Int64=100000, dt::Float64=.1)
	L = readdlm("data/ntw$(n)_lap_mat.csv",',')
	m = vec(readdlm("data/ntw$(n)_m.csv",','))
	d = vec(readdlm("data/ntw$(n)_d.csv",','))
	a = zeros(n)
	a[id] = a0
	f = zeros(n)
	f[id] = f0
	p = zeros(n)
	p[id] = p0

	Xs = generate_forced_time_series(L,m,d,(a,f,p),T,dt,.2*m)

	Df = 10
	n_period = 1.

	Lh,dh,ah,fh,ph = run_location_large_network(Xs, dt, Df, n_period)
	
	writedlm("data/ntw$(n)_afp_$id.csv",[ah fh ph],',')
end

tups = Array{Tuple{Int64,Int64,Float64,Float64,Float64,Int64,Float64},1}()
for n in [5,10,20]
	for id in 1:n
		push!(tups,(id,n,.2,.00945,pi/10,100000,.1))
	end
end

pmap(locate,tups)





