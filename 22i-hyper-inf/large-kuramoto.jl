using PyPlot, Distributed, Dates


n = 30; run = "003"

nkeep = 200
np = 2

addprocs(np)

@everywhere include("hyper_kuramoto.jl")
@everywhere include("hyper_inf.jl")

# #= Hypergraph from py ####################
A2 = zeros(n,n)
A2l = zeros(0,3)
A3 = zeros(n,n,n)
A3l = zeros(0,4)

el = readdlm("data/edgelist-n$n-"*run*".csv",',')
for l in 1:size(el)[1]
	global i,j,k,A2,A2l,A3,A3l
	i,j,k = el[l,:]
	if k == ""
		i,j = sort([i,j] .+ 1)
		A2[i,j] = A2[j,i] = 1.
		A2l = vcat(A2l,[i j 1.;j i 1.])
	else
		i,j,k = sort([i,j,k] .+ 1)
		A3[i,j,k] = A3[i,k,j] = A3[j,i,k] = A3[j,k,i] = A3[k,i,j] = A3[k,j,i] = 1.
		A3l = vcat(A3l,[i j k 1.;j i k 1.;k i j 1.])
	end
end
# =#


X = zeros(n,0)
Y = zeros(n,0)

δt = .01
maxiter = 20
for i in 1:100
	x,y = hyper_k(A2l,A3l,zeros(n),.1*rand(n),1.,0.,π/4,π/4,δt,maxiter,1e-6)
global 	X = [X x]
global 	Y = [Y y]
end
@info "Simulations done."

@info "============ STARTING INFERENCE ==========="
ooi = [2,3]
dmax = 2
coeff,ids = hyper_inf_par_filter(X,Y,ooi,dmax,nkeep,.01,.1)

A2this = zeros(0,3)
A3this = zeros(0,4)

d = get_d(n,dmax)

for i in 1:n
	c = vec(coeff[i])
	len = length(c)
	for l in collect(1:len)[abs.(c) .> 0.]
		id = ids[i][l]
		v0 = d[id,:]
		v = sort(setdiff(v0,[i,0]))
		if length(v) == 1 && 0 in v0
			global A2this = vcat(A2this,[i v[1] c[l]])
		elseif length(v) == 2
			global A3this = vcat(A3this,[i v[1] v[2] c[l]])
		end
	end
end

tpr2,fpr2,v,I0 = my_ROC_extended(abs.(A2this),A2l,n)
tpr3,fpr3,v,I0 = my_ROC_extended(abs.(A3this),A3l,n)
tpr,fpr,v,I0 = my_ROC_extended(abs.(A2this),A2l,abs.(A3this),A3l,n)

writedlm("data/rok-"*run*"-tpr2.csv",tpr2,',')
writedlm("data/rok-"*run*"-fpr2.csv",fpr2,',')
writedlm("data/rok-"*run*"-tpr3.csv",tpr3,',')
writedlm("data/rok-"*run*"-fpr3.csv",fpr3,',')
writedlm("data/rok-"*run*"-tpr.csv",tpr,',')
writedlm("data/rok-"*run*"-fpr.csv",fpr,',')


figure("ROK",(15,4))
subplot(1,3,1)
PyPlot.plot([0,1],[0,1],"--k",lw=.5)
PyPlot.plot(fpr,tpr)
xlabel("FPR")
ylabel("TPR")
title("ROC: all edges")
subplot(1,3,2)
PyPlot.plot([0,1],[0,1],"--k",lw=.5)
PyPlot.plot(fpr2,tpr2)
xlabel("FPR")
title("ROC: pairwise")
subplot(1,3,3)
PyPlot.plot([0,1],[0,1],"--k",lw=.5)
PyPlot.plot(fpr3,tpr3)
xlabel("FPR")
title("ROC: triadic")

PyPlot.savefig("rok-"*run*"-$(now()).pdf",format="pdf")



