using DelimitedFiles

include("centrality.jl")
include("kf.jl")
include("kuramoto_pert.jl")
include("res_dist.jl")

@info("Initialization")

Asp = readdlm("uk_adj_mat.csv",',') .+ 1.0
n = Int(maximum(Asp))
m = Int(size(Asp)[1]/2)

A = zeros(Int,n,n)
for i in 1:330
	A[Int(Asp[i,1]),Int(Asp[i,2])] = 1.0
end

line_list = Array{Array{Int,1},1}()
for i in 1:m
	push!(line_list,Array{Int,1}(vec(Asp[2*i-1,[1,2]])))
end

L = Array{Float64,2}(diagm(0 => vec(sum(A, dims = 2))) - A)

C = centralities(L)
Om = res_dist(L)

@info("Computing δKf1")

dKf1 = Array{Float64,1}()

#line_list = line_list[[1,20,40,60,80,100,120,140,160]]
#m = length(line_list)

for l in line_list
	b = -L[l[1],l[2]]
	Ome = Om[l[1],l[2]]
	if abs(b-Ome) < 1e-10
		push!(dKf1,Inf)
	else
		push!(dKf1,b*Ome/(1-b*Ome))
	end
end

ranked_line_idx = Array{Int,1}(vec(sortslices([dKf1 1:m],dims=1)[:,2]))
ranked_lines = line_list[ranked_line_idx]

#=
@info("Computing peformance measure(s)")

C1s = Array{Float64,1}()
C2s = Array{Float64,1}()

for l in line_list
	Ltemp = copy(L)
	Ltemp[l[1],l[2]] += 1
	Ltemp[l[2],l[1]] += 1
	Ltemp[l[1],l[1]] -= 1
	Ltemp[l[2],l[2]] -= 1
	
	C1,C2 = kuramoto_pert(Ltemp,zeros(n),zeros(n),("box",[.3,10,8]),1e-8,.1,10000)
	push!(C1s,C1)
	push!(C2s,C2)
end

C1_rank_idx = Array{Int,1}(vec(sortslices([C1s 1:m],dims=1)[:,2]))
C1_ranking = line_list[C1_rank_idx]
C2_rank_idx = Array{Int,1}(vec(sortslices([C2s 1:m],dims=1)[:,2]))
C2_ranking = line_list[C2_rank_idx]

x = vec(sortslices([ranked_line_idx 1:m],dims=1)[:,2])
y1 = vec(sortslices([C1_rank_idx 1:m],dims=1)[:,2])
figure()
PyPlot.plot(x,y1,"o")
xlabel("δKf1")
ylabel("C1")

y2 = vec(sortslices([C2_rank_idx 1:m],dims=1)[:,2])
figure()
PyPlot.plot(x,y2,"o")
xlabel("δKf1")
ylabel("C2")

figure()
PyPlot.plot(y1,y2,"o")
xlabel("C1")
ylabel("C2")


=#



