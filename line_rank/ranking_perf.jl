using DelimitedFiles,Statistics,LinearAlgebra,PyPlot,SparseArrays

include("kuramoto.jl")
include("L2B.jl")
include("uk_gen_idx.jl")
include("isconnected.jl")
include("rm_line.jl")
include("res_dist.jl")

lin_or_sin = "lin"
ntw = "uk"
plots = true
th_perf = true
if th_perf
	@warn "Theoretical P only works for constant ratio m/d"
end

taus = [0.025,0.075,0.125,0.175,0.225,0.275,0.325,0.375,0.425,0.475] # duration of the line contingency
taus = [0.6,]
P0 = .11

eps = 1e-6
max_iter = 1000000
h = 0.001

#=
if tau < h
	@warn("τ must be larger than h")
end
=#

# Choose network
if ntw == "uk"
	Asp = readdlm("uk_adj_mat.csv",',') 
elseif ntw == "ieee57"
	Asp = readdlm("ieee57_adj_mat.csv",',')
elseif ntw == "ieee57_unif"
	Asp = readdlm("ieee57_adj_mat.csv",',')
	Asp[:,3] = ones(size(Asp)[1])
end

# Initialize parameters
n = Int(maximum(Asp))
m = Int(size(Asp)[1]/2)

A = sparse(vec(Asp[:,1]),vec(Asp[:,2]),vec(Asp[:,3]))
L = spdiagm(0 => vec(sum(A,dims=2))) - A

line_list = Array{Tuple{Int64,Int64},1}()
for i in 1:m
	push!(line_list,(Int64(Asp[2*i-1,1]),Int64(Asp[2*i-1,2])))
end

mm = .2*(ones(n) + .6*rand(n) .- .3)
mm = .2*ones(n)
mm = zeros(n)
dd = .1*(ones(n) + .6*rand(n) .- .3)
dd = .1*ones(n)
M = spdiagm(0 => mm)
D = spdiagm(0 => dd)

P = zeros(n)
if ntw == "uk"
	P[gen_idx] = P0*ones(length(gen_idx))
elseif ntw == "ieee57"
	P = .002*vec(readdlm("P_57",','))
elseif ntw == "ieee57_unif"
	P = .002*vec(readdlm("P_57",','))
end
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

# Computing resistance distances
Om = res_dist(L)

@info("Computing δKf1")

dKf1 = Array{Float64,1}()

dist = Array{Float64,1}()
dist2 = Array{Float64,1}()
L2 = L^2
Om2 = res_dist(L2)
for l in line_list
	b = -L[l[1],l[2]]
	Ome = Om[l[1],l[2]]
	Ome2 = Om2[l[1],l[2]]
	if abs(b-Ome) < 1e-10
		push!(dKf1,Inf)
	else
		push!(dKf1,b*Ome/(1-b*Ome))
	end
	push!(dist,Ome)
	push!(dist2,Ome2)
end

ranked_line_idx = Array{Int,1}(vec(sortslices([dKf1 1:m],dims=1)[:,2]))
ranked_lines = line_list[ranked_line_idx]

Mi = 1
lM = 1
TM = 1
# for tau in taus
tau = taus[1]
	global L,M,D,Mi,lM,TM
	P1s = Array{Float64,1}()
	P2s = Array{Float64,1}()
	losses = Array{Float64,1}()
	P21s = Array{Float64,1}()
	P22s = Array{Float64,1}()
	Pths = Array{Float64,1}()
	co = 0
# Simulate line contingency for each line and compute performance measures
	for l in line_list
		global co
		co += 1
		@info("($(l[1]),$(l[2])), $co/$m")
#		global P1s,P2s,losses
		Ltemp = rm_line(L,l)
		if !isconnected(Ltemp)
			push!(P1s,Inf)
			push!(P2s,Inf)
			push!(P21s,Inf)
			push!(P22s,Inf)
			push!(losses,0.)
		else
			if lin_or_sin == "sin"
				xf,dxf,P1,P2,iter = kuramoto2_P(Ltemp,mm,dd,P,x1[1:n],x1[(n+1):(2*n)],x1[1:n],true,Int(ceil(tau/h)),eps,h)
				xff,dxff,P11,P22,iter = kuramoto2_P(L,mm,dd,P,xf[1:n],xf[(n+1):(2*n)],x1[1:n],true,max_iter,eps,h)
			elseif lin_or_sin == "lin"
				xf,dxf,P1,P2,iter = kuramoto2_lin_P(Ltemp,mm,dd,P,x1[1:n],x1[(n+1):(2*n)],x1[1:n],true,Int(ceil(tau/h)),eps,h)
				xff,dxff,P11,P22,iter = kuramoto2_lin_P(L,mm,dd,P,xf[1:n],xf[(n+1):(2*n)],x1[1:n],true,max_iter,eps,h)
			end
			
			push!(P1s,P1+P11)
			push!(P2s,P2+P22)
			push!(P21s,P2)
			push!(P22s,P22)
			if lin_or_sin == "sin"
				push!(losses,(L[l[1],l[2]]*sin(x1[l[1]]-x1[l[2]]))^2)
			elseif lin_or_sin == "lin"
				push!(losses,(L[l[1],l[2]]*(x1[l[1]]-x1[l[2]]))^2)
			end
		end

# Theoretical performance measure
		if th_perf
#= 
			X11 = zeros(n,n)
			gamma = dd[1]/mm[1]
			Mi = spdiagm(0 => 1 ./ mm)
			TM = eigvecs(Array(Mi.^(1/2)*Ltemp*Mi.^(1/2)))
			lM = diagm(0 => eigvals(Array(Mi.^(1/2)*Ltemp*Mi.^(1/2))))
			Gamma = sqrt.(Complex.(gamma^2 .- 4*lM))
			Mid = lM/(2*gamma) + exp(-gamma*tau)*((2*lM.^2)./(gamma*Gamma.^2) - (gamma*lM.*cosh.(tau*Gamma))./(2*Gamma.^2) - (lM.*sinh.(tau*Gamma))./(2*Gamma))
			
			X11 = TM*Mid*transpose(TM)
			
			B = x1[1:n]
			
			push!(Pths,transpose(B)*X11*B)
=#
			Omtemp = res_dist(Ltemp)
			push!(Pths,(2*tau^3)/(3*mm[1]*(1-L[l[1],l[2]]*Omtemp[l[1],l[2]])^2))
		end
	end
	
	
	
	
# Dealing with tree-like parts
	idx1 = Array{Int64,1}()
	idx2 = Array{Int64,1}()
	idx0 = Array{Int64,1}()
	for i in 1:m
		if dKf1[i] < 1e8
			push!(idx0,i)
		end
		if P1s[i]/losses[i] < 1e8
			push!(idx1,i)
		end
		if P2s[i]/losses[i] < 1e8
			push!(idx2,i)
		end
	end
	idx01 = intersect(idx0,idx1)
	idx02 = intersect(idx0,idx2)
	idx12 = intersect(idx1,idx2)
	idx = intersect(idx0,idx1,idx2)
	
# Compute correlations
	x = dKf1[idx0]
	y = P1s[idx1]./losses[idx1]
	z = P2s[idx2]./losses[idx2]
	t = dist
	
	x01 = dKf1[idx01]
	y01 = P1s[idx01]./losses[idx01]
	x02 = dKf1[idx02]
	z02 = P2s[idx02]./losses[idx02]
	y12 = P1s[idx12]./losses[idx12]
	z12 = P2s[idx12]./losses[idx12]
	t1 = dist[idx1]
	t2 = dist[idx2]
	
	r1 = sum((x01 .- mean(x01)).*(y01 .- mean(y01)))/(sqrt(sum((x01 .- mean(x01)).^2))*sqrt(sum((y01 .- mean(y01)).^2)))
	r2 = sum((x02 .- mean(x02)).*(z02 .- mean(z02)))/(sqrt(sum((x02 .- mean(x02)).^2))*sqrt(sum((z02 .- mean(z02)).^2)))
	r3 = sum((z12 .- mean(z12)).*(y12 .- mean(y12)))/(sqrt(sum((z12 .- mean(z12)).^2))*sqrt(sum((y12 .- mean(y12)).^2)))
	r4 = sum((t1 .- mean(t1)).*(y .- mean(y)))/(sqrt(sum((t1 .- mean(t1)).^2))*sqrt(sum((y .- mean(y)).^2)))
	r5 = sum((t2 .- mean(t2)).*(z .- mean(z)))/(sqrt(sum((t2 .- mean(t2)).^2))*sqrt(sum((z .- mean(z)).^2)))
	
# Plot
	if plots
		figure(ntw*"_"*lin_or_sin*"_contingency_$(P0)_$(tau)_$h",(15,8))
		
		PyPlot.subplot(1,5,1)
		PyPlot.plot(dKf1[idx01],P1s[idx01]./losses[idx01],"ob")
		xlabel("δKf1")
		ylabel("P1/losses")
		title("r = $(round(r1,digits=4))")
		
		PyPlot.subplot(1,5,2)
		PyPlot.plot(dKf1[idx02],P2s[idx02]./losses[idx02],"or")
		xlabel("δKf1")
		ylabel("P2/losses")
		title("r = $(round(r2,digits=4))")
		
		PyPlot.subplot(1,5,3)
		PyPlot.plot(P1s[idx12]./losses[idx12],P2s[idx12]./losses[idx12],"og")
		xlabel("P1/losses")
		ylabel("P2/losses")
		title("r = $(round(r3,digits=4))")
		
		PyPlot.subplot(1,5,4)
		PyPlot.plot(dist[idx1],P1s[idx1]./losses[idx1],"oy")
		xlabel("Ω")
		ylabel("P1/losses")
		title("r = $(round(r4,digits=4))")
		
		PyPlot.subplot(1,5,5)
		PyPlot.plot(dist[idx2],P2s[idx2]./losses[idx2],"om")
		xlabel("Ω")
		ylabel("P2/losses")
		title("r = $(round(r5,digits=4))")
		
		figure(ntw*"_"*lin_or_sin*"_contingency2_$(P0)_$(tau)_$h",(15,8))
		PyPlot.subplot(1,4,1)
		PyPlot.plot(P21s[idx2]./losses[idx2],"ok",label="Numerical P")
		PyPlot.plot(Pths[idx2]./losses[idx2],"xr",label="Theoretical P")
		title("P2 during cut")
		legend()
		PyPlot.subplot(1,4,2)
		PyPlot.plot(P22s[idx2]./losses[idx2],"ok")
		title("P2 after recovery")
		PyPlot.subplot(1,4,3)
		PyPlot.plot(dist2[idx1],P1s[idx1]./losses[idx1],"oy")
		title("P1 vs. Ω2")
		PyPlot.subplot(1,4,4)
		PyPlot.plot(dist2[idx2],P2s[idx2]./losses[idx2],"om")
		title("P2 vs. Ω2")
	end

	writedlm("data/P2s_$(P0)_$(tau).csv",P2s,',')
	writedlm("data/losses_$(P0)_$(tau).csv",losses,',')
#end



