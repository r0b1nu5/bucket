using Distributed

include("generate_time_series.jl")
include("final.jl")

#ntw = "ntw5"
#ntw = "ntw10"
ntw = "uk"

L = readdlm("data/"*ntw*"_lap_mat.csv",',')
m = vec(readdlm("data/"*ntw*"_m.csv",','))
Mi = diagm(0 => 1 ./ m)
d = vec(readdlm("data/"*ntw*"_d.csv",','))

Lm = Mi*L
dm = Mi*d

n = length(d)

a0 = .2
f0 = .009
p0 = pi/10

T = 100000
dt = .1

sig = mean(m)*ones(n)

tups = Array{Tuple{String,Array{Float64,2},Array{Float64,1},Array{Float64,1},Float64,Float64,Float64,Int64,Float64,Array{Float64,1},Int64},1}()

for i in 1:n
	push!(tups,(ntw,L,m,d,a0,f0,p0,T,dt,sig,i))
end

function plt(tup::Tuple{String,Array{Float64,2},Array{Float64,1},Array{Float64,1},Float64,Float64,Float64,Int64,Float64,Array{Float64,1},Int64})
	ntw,L,m,d,a0,f0,p0,T,dt,sig,i = tup

	n = length(dm)

	@info "$i/$n"

	a = zeros(n)
	a[i] = a0
	f = zeros(n)
	f[i] = f0
	p = zeros(n)
	p[i] = p0

	Xs = generate_forced_time_series(ntw, L, m, d, (a,f,p), T, dt, sig)
#	Xs = load_data(ntw, i)

	Lh,dh,ah,fh,ph = run_location_large_ntw(Xs, dt, 5)
#	Lh,dh,ah,fh,ph = run_location_small_ntw(Xs, dt)
	
	dL = Lh - Lm
	writedlm("data/"*ntw*"_reL_$(i).csv",sum(dL.^2)/sum(Lm.^2),',')
	
	da = ah[i] - a0
	writedlm("data/"*ntw*"_rea_$(i).csv",abs(da)/a0,',')
	
	am,id = findmax(ah)
	df = fh[id] - f0
	writedlm("data/"*ntw*"_ref_$(i).csv",abs(df)/f0,',')

	dp = ph[id] - p0
	writedlm("data/"*ntw*"_rep_$(i).csv",abs(dp)/p0,',')

	writedlm("data/"*ntw*"_ah_$(i).csv",ah,',')
end





