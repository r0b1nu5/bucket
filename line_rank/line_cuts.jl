using PyPlot,DelimitedFiles,SparseArrays,Statistics,LinearAlgebra

include("kuramoto.jl")
include("res_dist.jl")

ntw = "ieee118"
Lsp = readdlm(ntw*"_data/"*ntw*"_lap_mat.csv",',')
L = sparse(Lsp[:,1],Lsp[:,2],Lsp[:,3])

Om = res_dist(L)

n = size(L)[1]
m = round(Int,size(Lsp)[1]/2)

P0 = .001
P = P0*vec(readdlm(ntw*"_data/P_"*ntw*".csv",','))
P .-= mean(P)

g = .5
M = .2*ones(n)
D = g*M
LD = diagm(0 => D)^(-1/2)*L*diagm(0 => D)^(-1/2)

x0,dx0 = kuramoto2_lin(L,M,D,P,zeros(n),zeros(n))
th0 = x0[1:n]
om0 = x0[(n+1):(2*n)]

lines = Array{Tuple{Int64,Int64},1}()
nadirs = Array{Float64,2}(undef,n,0)
rocofs = Array{Float64,2}(undef,n,0)

dca = Array{Float64,2}(undef,4,0)

dth = Array{Float64,1}()
dua = Array{Float64,2}(undef,4,0)
rd = Array{Float64,1}()

full_measure = Array{Float64,2}(undef,4,0)

flow_m = Array{Float64,1}()

h = .01

for i in 1:n-1
	for j in i+1:n
		global L,n,M,D,P,th0,om0,h,nadirs,rocofs,lines,g,dca,dth,dua,rd,full_measure,flow_m
		@info "Line: ($i,$j)"
		if abs(L[i,j]) > 1e-6 && abs(1 + L[i,j]*Om[i,j]) > 1e-6
			Lr = copy(L)
			eij = zeros(n)
			eij[[i,j]] = [1,-1]
			Lr -= L[i,j]*eij*transpose(eij)
#			ls = eigvals(Array(Lr))
			
			Xs,dXs = kuramoto2_lin(Lr,M,D,P,th0,om0,true,false,Int(1e5),1e-6,h)
			xs = Xs[1:n,:]
			dxs = Xs[(n+1):(2*n),:]
			
			nadirs = [nadirs maximum(abs.(dxs),dims=2)]
			rocofs = [rocofs maximum(abs.(dxs[:,2:end]-dxs[:,1:end-1])./h,dims=2)]
			push!(lines,(i,j))
			
			LDr = diagm(0 => D)^(-1/2)*Lr*diagm(0 => D)^(-1/2)
			E = eigen(Array(LDr))
			lDs = E.values
			TD = E.vectors
			
			d2 = g^2-4*lDs[2]*g
			S2 = sqrt(abs(d2))
			if d2 >= 0
				F2 = (g+S2)/(g-S2)
				dc2 = (lDs[2]*g)/S2*F2^(-g/(2*S2))*(sqrt(F2)-sqrt(1/F2))
			else
				dc2 = (lDs[2]*g)/(S2)*exp(-(g*pi)/(2*S2))
			end
			
			d3 = g^2-4*lDs[3]*g
			S3 = sqrt(abs(d3))
			if d3 >= 0
				F3 = (g+S3)/(g-S3)
				dc3 = (lDs[3]*g)/S3*F3^(-g/(2*S3))*(sqrt(F3)-sqrt(1/F3))
			else
				dc3 = (lDs[3]*g)/(S3)*exp(-(g*pi)/(2*S3))
			end
			
			d4 = g^2-4*lDs[4]*g
			S4 = sqrt(abs(d4))
			if d4 >= 0
				F4 = (g+S4)/(g-S4)
				dc4 = (lDs[4]*g)/S4*F4^(-g/(2*S4))*(sqrt(F4)-sqrt(1/F4))
			else
				dc4 = (lDs[4]*g)/(S4)*exp(-(g*pi)/(2*S4))
			end
			
			dca = [dca [0.,dc2,dc3,dc4]]
			
			push!(dth,abs(th0[i]-th0[j]))
			dua = [dua vec(abs.(TD[i,1:4]-TD[j,1:4]))]
			push!(rd,Om[i,j])
			
			bij = -L[i,j]
			
			measure = (bij*dth[end])/(1-bij*Om[i,j])*(dua[:,end]./lDs[1:4])
			full_measure = [full_measure measure]
			
			push!(flow_m,abs(bij*dth[end])/min(M[i],M[j]))
		end
	end
end

figure("nadir")
PyPlot.plot(vec(maximum(nadirs,dims=1)),"xk")
PyPlot.plot(dca[2,:],"ob")
PyPlot.plot(dca[3,:],"or")
PyPlot.plot(dca[4,:],"og")
PyPlot.plot(full_measure[2,:],"sb")
PyPlot.plot(full_measure[3,:],"sr")
PyPlot.plot(full_measure[4,:],"sg")

figure("rocof")
PyPlot.plot(flow_m,"ob")
PyPlot.plot(vec(maximum(rocofs,dims=1)),"xr")








