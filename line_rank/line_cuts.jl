using PyPlot,DelimitedFiles,SparseArrays,Statistics,LinearAlgebra

include("kuramoto.jl")
include("res_dist.jl")
include("../default_color_cycle.jl")

h = .0001
col_num = 5

ntw = "ieee118"
Lsp = readdlm(ntw*"_data/"*ntw*"_lap_mat.csv",',')
L = sparse(Lsp[:,1],Lsp[:,2],Lsp[:,3])

Om = res_dist(L)

n = size(L)[1]
m = round(Int,size(Lsp)[1]/2)

P0 = .001
ratedP = vec(readdlm(ntw*"_data/P_"*ntw*".csv",','))
P = P0*ratedP
P .-= mean(P)

g = .5
# For inertias, see Laurent's Plos One, Appendix 2, Eq. (S2)
H = 2*pi .+ 3*rand(n) .- 1.5
omega0 = 50*2*pi
rP = 2*max.(abs.(ratedP),maximum(abs.(ratedP))/100*ones(n))
M = 2*H.*rP./omega0
D = g*M
LD = diagm(0 => D)^(-1/2)*L*diagm(0 => D)^(-1/2)

x0,dx0 = kuramoto2_lin(L,M,D,P,zeros(n),zeros(n))
th0 = x0[1:n]
om0 = x0[(n+1):(2*n)]

lines = Array{Tuple{Int64,Int64},1}()
nadirs = Array{Float64,2}(undef,n,0)
rocofs = Array{Float64,2}(undef,n,0)

dca = Array{Float64,2}(undef,n,0)

dth = Array{Float64,1}()
dua = Array{Float64,2}(undef,4,0)
rd = Array{Float64,1}()

full_measure = Array{Float64,2}(undef,4,0)

flow_m = Array{Float64,1}()

what = Array{Int64,1}()

idx = Array{Int64,1}()
scores = Array{Int64,2}(undef,n,0)

q1 = Array{Float64,1}()
q2 = Array{Float64,1}()
q3 = Array{Float64,1}()

c = 0

 #=
nadirs = readdlm("data/nadirs_"*ntw*"_$P0.csv",',')
rocofs = readdlm("data/rocofs_"*ntw*"_$P0.csv",',')
# =#

for i in 1:n-1
	for j in i+1:n
		global ntw,P0,L,n,M,D,P,th0,om0,h,nadirs,rocofs,lines,g,dca,dce,dth,dua,rd,full_measure,flow_m,c,what,idx,scores
		if abs(L[i,j]) > 1e-6 && abs(1 + L[i,j]*Om[i,j]) > 1e-6
			@info "Line: ($i,$j)"
			c += 1
			
			Lr = copy(L)
			eij = zeros(n)
			eij[[i,j]] = [1,-1]
			Lr -= L[i,j]*eij*transpose(eij)
#			ls = eigvals(Array(Lr))
			OmDr = res_dist(Lr)
			
			thf = pinv(Array(Lr))*P
 ##=			
			Xs,dXs = kuramoto2_lin(Lr,M,D,P,th0,om0,true,false,10,1e-6,h)
			xs = Xs[1:n,:]
			dxs = Xs[(n+1):(2*n),:]
			thf = xs[1:n,end]
			
			nadirs = [nadirs maximum(abs.(dxs),dims=2)]
			rocofs = [rocofs maximum(abs.(dxs[:,2:end]-dxs[:,1:end-1])./h,dims=2)]
			push!(lines,(i,j))
# =#
			
			LDr = diagm(0 => D)^(-1/2)*Lr*diagm(0 => D)^(-1/2)
			E = eigen(Array(LDr))
			tem = sortslices([E.values 1:n],dims=1)
			lDs = E.values[Int.(tem[:,2])]
			TD = E.vectors[:,Int.(tem[:,2])]
			Ss = Array{Float64,1}()
				
			temp = [0.,]
			Ss = [0.,]
			for k in 2:n
				d = g^2-4*lDs[k]*g
				S = sqrt(abs(d))
				push!(Ss,S)
				c0 = transpose(th0-thf)*sqrt(diagm(0 => D))*TD[:,k]
				if d >= 0
					F = (g+S)/(g-S)
					dc = (lDs[k]*g)/S*F^(-g/(2*S))*(sqrt(F)-sqrt(1/F))*c0
				else
					dc = (lDs[k]*g)/(S)*exp(-(g*pi)/(2*S))*c0
				end
				push!(temp,dc)
			end
			
			dca = [dca abs.(temp)]
	
			push!(dth,abs(th0[i]-th0[j]))
			dua = [dua vec(abs.(TD[i,1:4]-TD[j,1:4]))]
			push!(rd,Om[i,j])
			
			bij = -L[i,j]
			
			measure = (bij*dth[end])/(1-bij*Om[i,j])*(dua[:,end]./lDs[1:4])
			full_measure = [full_measure measure]
			
			push!(flow_m,abs(bij*dth[end])/min(M[i],M[j]))
			
			if maximum(dca[:,c]) > maximum(nadirs[:,c])
				@info "$c, $(maximum(dca[:,c])), $(maximum(nadirs[:,c]))"
				push!(what,c)
			end
			
			push!(idx,findmax(abs.(temp))[2])
			scores = [scores sortslices([(sortslices([abs.(temp) 1:n],dims=1))[:,2] 1:n],dims=1)]
			
			push!(q1,bij*abs(thf[i]-thf[j])/(1+bij*OmDr[i,j]))
			push!(q2,abs(TD[i,idx[end]]/D[i] - TD[j,idx[end]]/D[j]))
			push!(q3,Ss[idx[end]])
		end
	end
end


figure("rocof")
PyPlot.plot(flow_m,"ob")
PyPlot.plot(vec(maximum(rocofs,dims=1)),"xr")

figure("rocof 2")
PyPlot.plot(flow_m,vec(maximum(rocofs,dims=1)),"x",color=def_col[col_num])
PyPlot.plot([1e-8,10],[1e-8,10],"--k")
axis([0,maximum(flow_m),0,maximum(rocofs)])
xlabel("Initial flow on the removed line, normalized by inertia")
ylabel("Max. measured local RoCoF")

figure("rocof 3")
PyPlot.loglog(flow_m,vec(maximum(rocofs,dims=1)),"x",color=def_col[col_num])
PyPlot.plot([1e-8,10],[1e-8,10],"--k")
axis([minimum(flow_m),maximum(flow_m),minimum(flow_m),maximum(rocofs)])
xlabel("Initial flow on the removed line, normalized by inertia")
ylabel("Max. measured local RoCoF")

#=
figure("nadir")
for c in what
	PyPlot.plot([c-1,c-1],[0,dca[2,c]],"-r")
	global dca
end
for k in 1:n
	PyPlot.plot(dca[k,:],"o")
end
PyPlot.plot(vec(maximum(nadirs,dims=1)),"-xk",label="Max. measured nadir")
PyPlot.plot(vec(maximum(dca,dims=1)),"-xr",label="Max. \\dot{c}_a")
 #=
PyPlot.plot(full_measure[2,:],"sb")
PyPlot.plot(full_measure[3,:],"sr")
PyPlot.plot(full_measure[4,:],"sg")
=#
legend()
xlabel("Max. measured nadir")
ylabel("Max. \\dot{c}_a")

figure("nadir 2.1")
PyPlot.plot(vec(maximum(nadirs,dims=1)),vec(maximum(dca,dims=1)),"x",label="First \\dot{c}_a")
PyPlot.plot(vec(maximum(nadirs,dims=1)),[maximum(setdiff(dca[:,i],maximum(dca[:,i]))) for i in 1:size(dca)[2]],"x",label="Second \\dot{c}_a")
PyPlot.plot(vec(maximum(nadirs,dims=1)),[maximum(setdiff(dca[:,i],[maximum(setdiff(dca[:,i],maximum(dca[:,i]))),maximum(dca[:,i])])) for i in 1:size(dca)[2]],"x",label="Third \\dot{c}_a")
legend()
xlabel("Max. measured nadir")
ylabel("\\dot{c}_a")

figure("nadir 2.2")
PyPlot.loglog(vec(maximum(nadirs,dims=1)),vec(maximum(dca,dims=1)),"x",label="First \\dot{c}_a")
PyPlot.loglog(vec(maximum(nadirs,dims=1)),[maximum(setdiff(dca[:,i],maximum(dca[:,i]))) for i in 1:size(dca)[2]],"x",label="Second \\dot{c}_a")
PyPlot.loglog(vec(maximum(nadirs,dims=1)),[maximum(setdiff(dca[:,i],[maximum(setdiff(dca[:,i],maximum(dca[:,i]))),maximum(dca[:,i])])) for i in 1:size(dca)[2]],"x",label="Third \\dot{c}_a")
legend()
xlabel("Max. measured nadir")
ylabel("\\dot{c}_a")

figure("nadir 3")
PyPlot.bar((1:n).-.1,[sum(idx.==i) for i in 1:n]./length(idx),.2,label="Index has max. \\dot{c}_a")
PyPlot.bar((1:n).+.1,vec(sum(scores,dims=2))./sum(scores),.2,label="Score of each index")
legend()
xlabel("Eigenmode index")
ylabel("Proportion")

=#

##=
figure("nadir 4.1")
PyPlot.plot(q1,vec(maximum(nadirs,dims=1)),"x")
ylabel("Max. measured nadir")
xlabel("bij(thi-thj)/(1+bij*Omij)")

figure("nadir 4.2")
PyPlot.plot(q2,vec(maximum(nadirs,dims=1)),"x")
ylabel("Max. measured nadir")
xlabel("vi-vj")

figure("nadir 4.3")
PyPlot.plot(q3,vec(maximum(nadirs,dims=1)),"x")
ylabel("Max. measured nadir")
xlabel("Sa")

# =#


