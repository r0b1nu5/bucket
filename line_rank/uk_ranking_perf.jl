using DelimitedFiles,Statistics,LinearAlgebra,PyPlot

include("kuramoto.jl")
include("L2B.jl")
include("uk_gen_idx.jl")
include("isconnected.jl")
include("rm_line.jl")
include("res_dist.jl")

lin_or_sin = "sin"
ntw = "uk"

tau = 0.25 # duration of the line contingency
P0 = .2

eps = 1e-6
max_iter = 100000
h = 0.025

if tau < h
	@warn("τ must be larger than h")
end

# Choose network
if ntw == "uk"
	Asp = readdlm("uk_adj_mat.csv",',') .+ 1
end

# Initialize parameters
n = Int(maximum(Asp))
m = Int(size(Asp)[1]/2)

A = sparse(vec(Asp[:,1]),vec(Asp[:,2]),ones(size(Asp)[1]))
L = spdiagm(0 => vec(sum(A,dims=2))) - A

line_list = Array{Tuple{Int64,Int64},1}()
for i in 1:m
	push!(line_list,(Int64(Asp[2*i-1,1]),Int64(Asp[2*i-1,2])))
end

mm = .2*ones(n)
dd = .1*ones(n)

P = zeros(n)
P[gen_idx] = P0*ones(length(gen_idx))
P .-= mean(P)

th0 = zeros(n)
omeg0 = zeros(n)
x0 = [th0;omeg0]

# Find fixed point
if lin_or_sin == "sin"
	xs1,dxs1 = kuramoto2(L,mm,dd,P,x0[1:n],x0[(n+1):(2*n)])
elseif lin_or_sin == "lin"
	xs1,dxs1 = kuramoto2_lin(L,mm,dd,P,x0[1:n],x0[(n+1):(2*n)])
end
x1 = vec(xs1[:,end])


P1s = Array{Float64,1}()
P2s = Array{Float64,1}()

# Simulate line contingency for each line and compute performance measures
for l in line_list
	@info(l)
	global P1s,P2s
	Ltemp = rm_line(L,l)
	if !isconnected(Ltemp)
		push!(P1s,Inf)
		push!(P2s,Inf)
	else
		if lin_or_sin == "sin"
			xs2,dxs2 = kuramoto2(Ltemp,mm,dd,P,x1[1:n],x1[(n+1):(2*n)],true,true,Int(ceil(tau/h)))
			xs3,dxs3 = kuramoto2(L,mm,dd,P,vec(xs2[1:n,end]),vec(xs2[(n+1):(2*n),end]),true)
		elseif lin_or_sin == "lin"
			xs2,dxs2 = kuramoto2_lin(Ltemp,mm,dd,P,x1[1:n],x1[(n+1):(2*n)],true,true,Int(ceil(tau/h)))
			xs3,dxs3 = kuramoto2_lin(L,mm,dd,P,vec(xs2[1:n,end]),vec(xs2[(n+1):(2*n),end]),true)
		end
		
		ths = [xs2[1:n,:] xs3[1:n,:]]
		omegs = [xs2[(n+1):(2*n),:] xs3[(n+1):(2*n),:]]
		
		if lin_or_sin == "sin"		
			push!(P1s,h*sum((ths-repeat(x1[1:n],outer=(1,size(ths)[2]))).^2)/abs(1-cos(x1[l[1]]-x1[l[2]])))
			push!(P2s,h*sum(omegs.^2)/abs(1-cos(x1[l[1]]-x1[l[2]])))
		elseif lin_or_sin == "lin"
			push!(P1s,h*sum((ths-repeat(x1[1:n],outer=(1,size(ths)[2]))).^2)/(x1[l[1]]-x1[l[2]])^2)
			push!(P2s,h*sum(omegs.^2)/(x1[l[1]]-x1[l[2]])^2)
		end
			
	end
end

# Computing resistance distances
Om = res_dist(L)

line_list = Array{Array{Int,1},1}()
dist = Array{Float64,1}()
for i in 1:m
	push!(line_list,Array{Int,1}(vec(Asp[2*i-1,[1,2]])))
	push!(dist,Om[line_list[end][1],line_list[end][2]])
end

@info("Computing δKf1")

dKf1 = Array{Float64,1}()

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


# Dealing with tree-like parts
idx = Array{Int64,1}()
for i in 1:m
	if !(isinf(dKf1[i]) || isinf(P1s[i]) || isinf(P2s[i]))
		push!(idx,i)
	end
end

# Compute correlations
x = dKf1[idx]
y = P1s[idx]
z = P2s[idx]
t = dist[idx]

r1 = sum((x .- mean(x)).*(y .- mean(y)))/(sqrt(sum((x .- mean(x)).^2))*sqrt(sum((y .- mean(y)).^2)))
r2 = sum((x .- mean(x)).*(z .- mean(z)))/(sqrt(sum((x .- mean(x)).^2))*sqrt(sum((z .- mean(z)).^2)))
r3 = sum((z .- mean(z)).*(y .- mean(y)))/(sqrt(sum((z .- mean(z)).^2))*sqrt(sum((y .- mean(y)).^2)))
r4 = sum((t .- mean(t)).*(y .- mean(y)))/(sqrt(sum((t .- mean(t)).^2))*sqrt(sum((y .- mean(y)).^2)))
r5 = sum((t .- mean(t)).*(z .- mean(z)))/(sqrt(sum((t .- mean(t)).^2))*sqrt(sum((z .- mean(z)).^2)))

# Plot
figure(ntw*"_"*lin_or_sin*"_contingency_$P0",(15,8))

PyPlot.subplot(1,5,1)
PyPlot.plot(dKf1,P1s,"ob")
xlabel("δKf1")
ylabel("P1")
title("r = $(round(r1,digits=4))")

PyPlot.subplot(1,5,2)
PyPlot.plot(dKf1,P2s,"or")
xlabel("δKf1")
ylabel("P2")
title("r = $(round(r2,digits=4))")

PyPlot.subplot(1,5,3)
PyPlot.plot(P1s,P2s,"og")
xlabel("P1")
ylabel("P2")
title("r = $(round(r3,digits=4))")

PyPlot.subplot(1,5,4)
PyPlot.plot(dist,P1s,"oy")
xlabel("Ω")
ylabel("P1")
title("r = $(round(r4,digits=4))")

PyPlot.subplot(1,5,5)
PyPlot.plot(dist,P2s,"om")
xlabel("Ω")
ylabel("P2")
title("r = $(round(r5,digits=4))")



