using PyPlot,SparseArrays,Statistics,DelimitedFiles,LinearAlgebra,Distributions

include("kuramoto.jl")
include("res_dist.jl")
include("rm_line.jl")
include("L2B.jl")
include("isconnected.jl")

thr = 2206
h = 1e-4


n_simu = 1000

ntw = "ieee118"
Lsp = readdlm(ntw*"_data/"*ntw*"_lap_mat.csv",',')
L = sparse(Lsp[:,1],Lsp[:,2],Lsp[:,3])
Li = pinv(Array(L))

Bsp = readdlm(ntw*"_data/"*ntw*"_inc_mat.csv",',')

n = size(L)[1]
m = round(Int,(size(Lsp)[1]-n)/2)

P0 = .01
ratedP = vec(readdlm(ntw*"_data/P_"*ntw*".csv",','))
P = P0*ratedP
P .-= mean(P)

g = .5
# #=
H = 2*pi .+ 2*rand(n) .- 1.
omega0 = 50*2*pi
rP = 2*max.(abs.(ratedP),maximum(abs.(ratedP))/100*ones(n))
M = 2*H.*rP./omega0

H = 2*pi .+ 2*rand(n) .- 1.
omega0 = 50*2*pi
rP = 2*max.(abs.(ratedP),maximum(abs.(ratedP))/100*ones(n))
D = H.*rP./omega0
# =#
#M = vec(readdlm(ntw*"_data/M.csv",','))
#D = g*M

cutable_lines = Array{Tuple{Int64,Int64},1}()
for i in 1:n-1
	for j in i+1:n
		if abs(L[i,j]) > 1e-4 && isconnected(rm_line(L,(i,j)))
			if M[i] < M[j]
				push!(cutable_lines,(i,j))
			else
				push!(cutable_lines,(j,i))
			end
		end
	end
end

E = eigen(Array(L))
ls = E.values
us = E.vectors

varP = abs.(P)/3

Ps = zeros(n,n_simu)
for i in 1:n
	Ptemp = rand(Normal(P[i],varP[i]),1,n_simu)
	Ps[i,:] = Ptemp
end
for i in 1:n_simu
	Ps[:,i] .-= mean(Ps[:,i])
end
EP = sum(Ps,dims=2)./n_simu
covP = cov(Ps,dims=2)
#Ps = Array{Float64,2}(undef,n,0)

Es = Array{Float64,2}(undef,3,0)
vars = Array{Float64,2}(undef,3,0)

c = 0

for l in cutable_lines
	global c,L,Li,P,M,n,us,Ps,ls,varP,n_simu,h,Es,vars,covP
	c += 1
	@info "Line $c/$(length(cutable_lines))"
	
	i = l[1]
	j = l[2]
	eij = zeros(n)
	eij[[i,j]] = [1,-1]
#=
	P2use = Array{Int64,1}()
	for sim in 1:n_simu
		if (transpose(eij)*Li*Ps[:,simu])[1] > 0
			push!(P2use,sim)
		end
	end
	
	Px = Ps[:,P2use]
	n2use = length(P2use)
	EP = sum(Px,dims=2)./n2use
=#
	
	th = Li*EP
	
#	Ei_p = -L[i,j]/M[i]*abs(th[i]-th[j])
#	Ej_p = M[i]/M[j]*Ei_p

	Ei_p = -L[i,j]/M[i]*(th[i]-th[j])
	Ej_p = -M[i]/M[j]*Ei_p
		
	vari_p = 0.
	
	for a in 2:n
		for b in 2:n
	#		vari_p += (L[i,j])^2/M[i]^2*(us[i,a]-us[j,a])*(us[i,b]-us[j,b])/(ls[a]*ls[b])*(transpose(us[:,a])*diagm(0 => varP)*us[:,b])[1]
			vari_p += (L[i,j])^2/M[i]^2*(us[i,a]-us[j,a])*(us[i,b]-us[j,b])/(ls[a]*ls[b])*(transpose(us[:,a])*covP*us[:,b])[1]
		end
	end
	
	varj_p = M[i]^2/M[j]^2*vari_p
	
	Lr = rm_line(L,l)
	
	rcfs = Array{Float64,1}()
	for sim in 1:n_simu
		@info "Simu $sim"
#=	
		dP = Array{Float64,1}()
		for i in 1:n
			push!(dP,rand(Normal(0,abs(P[i])/3)))
		end
		dP .-= mean(dP)
		Ps = [Ps P+dP]
=#		
		th0 = Li*vec(Ps[:,sim])
		
		rcf,iidd = rocof_kuramoto2_lin(Lr,M,D,vec(Ps[:,sim]),th0,zeros(n),false,10,1e-4,h)
		
		push!(rcfs,rcf)
		if iidd > 1
			@warn "Max RoCoF at iter = $iidd"
		end
	end

	E_m = mean(rcfs)
	
	var_m = var(rcfs)
	
	Es = [Es [Ei_p,Ej_p,E_m]]
	vars = [vars [vari_p,varj_p,var_m]]
end
	
 #= 	
figure()
PyPlot.plot(Es[1,:],"o",label="Node i")
PyPlot.plot(Es[2,:],"o",label="Node j")
PyPlot.plot(Es[3,:],"x",label="Max. measured")
title("Mean RoCoFs")
legend()

figure()
PyPlot.plot(Es[1,:],Es[3,:],"x",label="Node i")
PyPlot.plot(Es[2,:],Es[3,:],"x",label="Node j")
xlabel("Predicted expectation")
ylabel("Measured expectation")
legend()

figure()
PyPlot.plot(vars[1,:],"o",label="Node i")
PyPlot.plot(vars[2,:],"o",label="Node j")
PyPlot.plot(vars[3,:],"x",label="Max. measured")
title("RoCoFs' variances")
legend()

figure()
PyPlot.plot(vars[1,:],vars[3,:],"x",label="Node i")
PyPlot.plot(vars[2,:],vars[3,:],"x",label="Node j")
xlabel("Predicted variance")
ylabel("Measured variance")
legend()
# =#

 ##=
writedlm("data/Es_$(thr).csv",Es,',')
writedlm("data/vars_$(thr).csv",vars,',')
# =#









