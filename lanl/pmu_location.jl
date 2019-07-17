using DelimitedFiles, LightGraphs, SimpleWeightedGraphs, Dates, Random

include("locate_forced.jl")
include("res_dist.jl")
include("load_uk.jl")

n_sample = 100
n_pmu = 20

T = 50000
dt = .1
f = 2e-4
w = 2*pi*f
forced_id = 45

g = SimpleWeightedGraph(Array{Int64,1}(Asp[:,1]),Array{Int64,1}(Asp[:,2]),Asp[:,3])
gdist = zeros(n,n)
for i in 1:n-1
	dsp = dijkstra_shortest_paths(g,i)
	for j in i+1:n
		gdist[i,j] = length(enumerate_paths(dsp,j)) - 1
		gdist[j,i] = gdist[i,j]
	end
end

rdist = res_dist(L)

Xs = readdlm("data/uk45_2_forced_inertialess_0.0002_50000_0.1.csv",',')

means_gd = Array{Float64,1}()
means_rd = Array{Float64,1}()
ranks_full = Array{Float64,1}()
ranks_slow = Array{Float64,1}()
idxs = Array{Int64,2}(undef,n_pmu,0)

sample_num = 0

while sample_num < n_sample
	global sample_num, n_pmu, n, gdist, rdist, Xs, L, w, T, dt, forced_id, means_gd, means_rd, ranks_full, ranks_slow, idxs
	
	sample_num += 1
	@info "$(now()) -- Sample $(sample_num)"
	
	pmu_idx = shuffle(1:n)[1:n_pmu]
	
	mean_gd = sum(gdist[pmu_idx,pmu_idx])/(n_pmu^2-n_pmu)
	mean_rd = sum(rdist[pmu_idx,pmu_idx])/(n_pmu^2-n_pmu)
		
	X = Xs[pmu_idx,:]
	
	err_id_full = locate_forcing_full_gen(X,pmu_idx,L,w,T,dt)
	err_id_slow = locate_forcing_slow(X,pmu_idx,L)
	
	rank_full = sortslices([abs.(err_id_full[:,2] .- forced_id) 1:n],dims=1)[1,2]
	rank_slow = sortslices([abs.(err_id_slow[:,2] .- forced_id) 1:n],dims=1)[1,2]
	
	push!(means_gd,mean_gd)
	push!(means_rd,mean_rd)
	push!(ranks_full,rank_full)
	push!(ranks_slow,rank_slow)
	idxs = [idxs pmu_idx]
end

writedlm("data/mgd.csv",means_gd,',')
writedlm("data/mrd.csv",means_rd,',')
writedlm("data/rf.csv",ranks_full,',')
writedlm("data/rs.csv",ranks_slow,',')
writedlm("data/pmu_idxs.csv",idxs,',')

figure("dist vs. rank")

subplot(221)
PyPlot.plot(means_gd,ranks_full,"x")
xlabel("mean geo. dist. btw PMUs")
ylabel("rank attributed to the forced node (full_gen code)")

subplot(223)
PyPlot.plot(means_gd,ranks_slow,"x")
xlabel("mean geo. dist. btw PMUs")
ylabel("rank attributed to the forced node (slow_code)")

subplot(222)
PyPlot.plot(means_rd,ranks_full,"x")
xlabel("mean res. dist. btw PMUs")
ylabel("rank attributed to the forced node (full_gen code)")

subplot(224)
PyPlot.plot(means_rd,ranks_slow,"x")
xlabel("mean res. dist. btw PMUs")
ylabel("rank attributed to the forced node (slow code")







