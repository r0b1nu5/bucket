using DelimitedFiles,Statistics,PyPlot,SparseArrays,LinearAlgebra

include("kuramoto.jl")
include("res_dist.jl")
include("rm_line.jl")
include("L2B.jl")

n_simu = 100

ntw = "ieee118"
Lsp = readdlm(ntw*"_data/"*ntw*"_lap_mat.csv",',')
L = sparse(Lsp[:,1],Lsp[:,2],Lsp[:,3])
Li = pinv(Array(L))

B,w = L2B(L)

Om = res_dist(L)

n = size(L)[1]
m = round(Int,(size(Lsp)[1]-n)/2)

P0 = .001
ratedP = vec(readdlm(ntw*"_data/P_"*ntw*".csv",','))
P = P0*ratedP
P .-= mean(P)

g = .5
# For inertias, see Laurent's Plos One, Appendix 2, Eq. (S2)
H = 2*pi .+ 2*rand(n) .- 1.
omega0 = 50*2*pi
rP = 2*max.(abs.(ratedP),maximum(abs.(ratedP))/100*ones(n))
M = 2*H.*rP./omega0
D = g*M
LD = diagm(0 => D)^(-1/2)*L*diagm(0 => D)^(-1/2)

E = eigen(Array(L))
TD = E.vectors
ls = E.values

X = transpose(B)*TD./repeat(transpose(ls),m,1)





# Lines to cut: largest bij/mi, median bij/mi, smallest bij/mi,...

#lines2cut = [(5,4),(37,34),(12,16),(96,94),(89,85),(69,49)]
lines2cut = [(65,66),(68,81),(53,54),(64,61)]
# Distribution P: large variance on slow modes, large variance on fast modes,...

id2 = [100,92,88,89,101,102,90,91,86,103,87,104,106,105,107,108,109,110,112,111]
id3 = Array{Int,1}()#[60,109,47,49,110,45,59,48,46,112,111,50,54,55,56,57,51,58,53,52]
id4 = Array{Int,1}()#[92,103,83,104,91,89,106,90,105,107,88,84,85,108,109,110,112,111,86,87]
id234 = union(id2,id3,id4)
n234 = length(id234)

id118 = [78,67,60,79,37,47,70,75,77,61,63,38,80,66,64,69,81,65,116,68]
id117 = Array{Int,1}()#[2,113,34,26,37,13,38,17,1,10,12,7,30,11,9,3,6,8,4,5]
id116 = Array{Int,1}()#[4,8,5,66,41,17,116,64,65,30,19,43,40,33,39,38,35,36,34,37]
id876 = union(id118,id117,id116)
n876 = length(id876)

th0 = Li*P


sim = 0
h = .01

rcfs1 = Array{Float64,2}(undef,length(lines2cut),0)
rcfs3 = Array{Float64,2}(undef,length(lines2cut),0)

while sim < n_simu
	global sim,n,n234,n876,id234,id876,P,lines2cut,L,M,D,th0,h,rcfs1,rcfs3
	sim += 1
	
	dp1 = (2*rand(n234) .- 1).*.1
	dp1 .-= mean(dp1)
	
	dP1 = zeros(n)
	dP1[id234] = dp1
	
	th1 = Li*(P+dP1)
	
	dp2 = (2*rand(n876) .- 1).*.1
	dp2 .-= mean(dp2)
	
	dP2 = zeros(n)
	dP2[id876] = dp2
	
	th2 = Li*(P+dP2)
	
	rcf1 = Array{Float64,1}()
	rcf3 = Array{Float64,1}()
	
	for l in lines2cut
		Lr = rm_line(L,l)
		@info "Simu $sim, Line $(l)"
		
		rcf = rocof_kuramoto2_lin(Lr,M,D,P+dP1,th1,zeros(n),false,Int(1e5),1e-4,h)
#=
		xs = Xs[1:n,:]
		dxs = Xs[(n+1):(2*n),:]
		
		r,i = findmax(abs.(dxs[l[1],2:end]-dxs[l[1],1:end-1])./h)
		push!(rcf1,r)
		if i > 1
			@info "1) Largest RoCoF at time = $i"
		end
#		push!(rcf1,maximum(abs.(dxs[l[1],2:end]-dxs[l[1],1:end-1])./h))		
#		rcf2 = maximum(abs.(dxs[l[2],2:end]-dxs[l[2],1:end-1])./h)
=#
		push!(rcf1,rcf)
		
		rcf = rocof_kuramoto2_lin(Lr,M,D,P+dP2,th2,zeros(n),false,Int(1e5),1e-4,h)
#= 
		xs = Xs[1:n,:]
		dxs = Xs[(n+1):(2*n),:]
		
		r,i = findmax(abs.(dxs[l[1],2:end]-dxs[l[1],1:end-1])./h)
		push!(rcf3,r)
		if i > 1
			@info "2) Largest RoCoF at time = $i"
		end	
#		push!(rcf3,maximum(abs.(dxs[l[1],2:end]-dxs[l[1],1:end-1])./h))
#		rcf4 = maximum(abs.(dxs[l[2],2:end]-dxs[l[1],1:end-1])./h)
=#
		push!(rcf3,rcf)
	end
	
	rcfs1 = [rcfs1 rcf1]
	rcfs3 = [rcfs3 rcf3]
end

c = 0

for l in lines2cut
	global c,rcfs1,rcfs3,n_simu
	c += 1
	
	m1 = round(mean(rcfs1[c,:]),digits=6)
	sd1 = round(sqrt(var(rcfs1[c,:])),digits=6)
	m3 = round(mean(rcfs3[c,:]),digits=6)
	sd3 = round(sqrt(var(rcfs3[c,:])),digits=6)
	
	figure("$(l)")
	
	subplot(121)
	PyPlot.hist(rcfs1[c,:],ceil(Int,n_simu/10),color="blue")
	PyPlot.plot(m1,.02*n_simu,"ok",label="μ=$m1")
	PyPlot.plot([m1-sd1,m1+sd1],[.02*n_simu,.02*n_simu],"--k",label="σ=$sd1")
	xlabel("RoCoF (large variance on slow modes)")
	legend()
	
	subplot(122)
	PyPlot.hist(rcfs3[c,:],ceil(Int,n_simu/10),color="green")
	PyPlot.plot(m3,.02*n_simu,"ok",label="μ=$m3")
	PyPlot.plot([m3-sd3,m3+sd3],[.02*n_simu,.02*n_simu],"--k",label="σ=$sd3")
	xlabel("RoCoF (large variance on fast modes)")
	legend()
end
	
	
	
	
	
	
		






