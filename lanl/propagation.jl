using DelimitedFiles, PyPlot, FFTW, LightGraphs, SimpleWeightedGraphs

include("generate_time_series.jl")
include("res_dist.jl")

fs = round.(exp.(LinRange(log(.01),log(8),20)),digits=2)
idx = [112, 65, 76] # 112 is on u_2\u_n, 65 is on u_n, 76 is neither on u_2 or u_n
n = 120

Lsp = readdlm("data/lap_mat.csv",',')
Asp = readdlm("data/uk_adj_mat.csv",',')
L = zeros(n,n)
for i in 1:size(Lsp)[1]
	L[Int(Lsp[i,1]),Int(Lsp[i,2])] = Lsp[i,3]
end
g = SimpleWeightedGraph(Array{Int64,1}(Asp[:,1]),Array{Int64,1}(Asp[:,2]),Asp[:,3])

Om = res_dist(L)

m = ones(n)
d = ones(n)

dt = .05
T = 20000

sig = ones(n)

generate = false

if generate
	for id in idx
		for fr in fs
			global L,m,d,T,dt,sig
			c = zeros(n)
			c[id] = 1.
			f = zeros(n)
			f[id] = fr
			phi = zeros(n)
			forc = (c,f,phi)
			generate_forced_time_series("uk$(id)",L,m,d,forc,T,dt,sig)
		end
	end
end

# #=
d = Array{Int64,2}(undef,length(fs),0)
rd = Array{Float64,2}(undef,length(fs),0)
for id in idx
	global d,rd,Om
	dist_i = Array{Int,1}()
	mnodes = Array{Int,1}()
	rdist_i = Array{Float64,1}()
	rmnodes = Array{Int,1}()
	for fr in fs
		@info "i = $(id), f = $(fr)"
		Xs = readdlm("data/uk$(id)_forced_$(fr)_$(T)_$(dt).csv",',')
		
		k = round(Int,fr*dt*T) + 1
		
		reached = [id,]
		dist_ij = [0,]
		rdist_ij = [0.,]
		
		for i in 1:n
			fX = abs.(real(fft(Xs[i,:])))
			nx = fX[[k.-(8:-1:1);k.+(1:8)]]
			dx = maximum(nx) - minimum(nx)
			if fX[k] > maximum(nx) + dx
				push!(reached,i)
				push!(dist_ij,length(enumerate_paths(dijkstra_shortest_paths(g,id),i)) - 1)
				push!(rdist_ij,Om[id,i])
			end
		end
		mad,maid = findmax(dist_ij)
		push!(mnodes,reached[maid])
		push!(dist_i,mad)
		rmad,rmaid = findmax(rdist_ij)
		push!(rmnodes,reached[rmaid])
		push!(rdist_i,maximum(rmad))
	end
	d = [d dist_i]
	rd = [rd rdist_i]
	
	figure("propagation 1")
	PyPlot.plot(fs,dist_i,"-o",label="idx = $(id)")
	for i in 1:length(fs)
		PyPlot.text(fs[i],dist_i[i],"$(mnodes[i])")
	end
	figure("propagation 2")
	PyPlot.plot(fs,rdist_i,"-o",label="idx = $(id)")
	for i in 1:length(fs)
		PyPlot.text(fs[i],rdist_i[i],"$(rmnodes[i])")
	end

end

figure("propagation 1")
legend()
xlabel("frequency")
ylabel("distance")	

figure("propagation 2")
legend()
xlabel("frequency")
ylabel("distance")	
			





# =#





