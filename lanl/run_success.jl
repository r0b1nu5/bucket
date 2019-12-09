# Runs the algorithm n_run times for a forcing located at each node of the network.

using Distributed

include("generate_time_series.jl")
include("final.jl")

go = false
if go
	@info "Computation will run."
else
	@info "Computation will not run. To run it set 'go' to 'true'"
end

n_run = 100

ntw = "uk10"
L = readdlm("data/"*ntw*"_lap_mat.csv",',')
m = vec(readdlm("data/"*ntw*"_m.csv",','))
d = vec(readdlm("data/"*ntw*"_d.csv",','))
	
n = length(d)
	
a0 = .2
f0 = .001
p0 = pi/10
	
T = 100000
dt = .1
	
sig = mean(m)*ones(n)

tups = Array{Tuple{String,Array{Float64,2},Array{Float64,1},Array{Float64,1},Float64,Float64,Float64,Int64,Float64,Array{Float64,1},Int64,Int64},1}()

for i in 1:n_run
	for j in 1:n
		push!(tups,(ntw,L,m,d,a0,f0,p0,T,dt,sig,j,i))
	end
end

function plts_success_rate(tup::Tuple{String,Array{Float64,2},Array{Float64,1},Array{Float64,1},Float64,Float64,Float64,Int64,Float64,Array{Float64,1},Int64,Int64})
	ntw,L,m,d,a0,f0,p0,T,dt,sig,node,run = tup
	
	n = length(m)

	@info "Run $run, node $node/$n"

	a = zeros(n)
	a[node] = a0
	f = zeros(n)
	f[node] = f0
	p = zeros(n)
	p[node] = p0

	Mi = diagm(0 => 1 ./ m)
	Lm = Mi*L
	
	Xs = generate_forced_time_series(ntw*"_node$(node)_run$(run)", L, m, d, (a,f,p), T, dt, sig)
#	Xs = load_data(ntw, i)
	
	if n < 10
		Lh,dh,ah,fh,ph = run_location_small_ntw(Xs,dt)
	else
		Lh,dh,ah,fh,ph = run_location_large_ntw(Xs, dt, 5)
	end
	
#=
	dL = Lh - Lm
	writedlm("data/"*ntw*"_reL_$(i).csv",sum(dL.^2)/sum(Lm.^2),',')
	
	da = ah[i] - a0
	writedlm("data/"*ntw*"_rea_$(i).csv",abs(da)/a0,',')
	
	am,id = findmax(ah)
	df = fh[id] - f0
	writedlm("data/"*ntw*"_ref_$(i).csv",abs(df)/f0,',')

	dp = ph[id] - p0
	writedlm("data/"*ntw*"_rep_$(i).csv",abs(dp)/p0,',')
=#

	writedlm("data/"*ntw*"_$(node)_$(run)_ah.csv",ah,',')
end


if go
	pmap(plts_success_rate,tups)
end


