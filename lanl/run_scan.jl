# Runs the algorithm n_run times for a forcing located at node id for each value of the parameter 'param' in the range LinRange(pmin,pmax,np).

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
id = 76

param = "fre" # ["amp","fre","pha","tim","ste","noi"]
pmin = .0005
pmax = .007
np = 14

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

tups = Array{Tuple{String,String,Float64,Float64,Int64,Int64,Array{Float64,2},Array{Float64,1},Array{Float64,1},Float64,Float64,Float64,Int64,Float64,Array{Float64,1},Int64,Int64},1}()

for i in 1:n_run
	for j in 1:np
		push!(tups,(ntw,param,pmin,pmax,np,j,L,m,d,id,i))
	end
end

function plts_scan(tup::Tuple{String,String,Float64,Float64,Int64,Int64,Array{Float64,2},Array{Float64,1},Array{Float64,1},Int64,Int64})
	ntw,param,pmin,pmax,np,val,L,m,d,a0,f0,p0,T,dt,sig,node,run = tup
	
	n = length(m)

	@info "Run $run for $param, value $j/$np"

	a = zeros(n)
	a[node] = .2
	f = zeros(n)
	f[node] = .001
	p = zeros(n)
	p[node] = pi/10
	T = 100000
	dt = .1
	sig = mean(m)*ones(n)

	if param == "amp"
		a[node] = LinRange(pmin,pmax,np)[val]
	elseif param == "fre"
		f[node] = LinRange(pmin,pmax,np)[val]
	elseif param == "pha"
		p[node] = LinRange(pmin,pmax,np)[val]
	elseif param == "tim"
		T = round(Int,LinRange(pmin,pmax,np)[val])
	elseif param == "ste"
		dt = LinRange(pmin,pmax,np)[val]
	elseif param == "noi"
		sig = ones(n)*LinRange(pmin,pmax,np)[val]
	end
	
	Mi = diagm(0 => 1 ./ m)
	Lm = Mi*L
	
	Xs = generate_forced_time_series(ntw*"_$(param)($pmin,$pmax,$np,$val)_run$(run)", L, m, d, (a,f,p), T, dt, sig)
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

	writedlm("data/"*ntw*"_"*param*"($pmin,$pmax,$np,$val)_run$(run)_ah.csv",ah,',')
end


if go
	pmap(plts_success_rate,tups)
end


