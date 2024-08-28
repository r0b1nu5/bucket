using PyPlot

include("hyper_kuramoto.jl")
include("hyper_inf.jl")
include("gen_rand_hyperg.jl")

n = 300

 #= ################ ER Hypergraph ######################
p2 = .4
m2 = round(Int64,n*(n-1)/2*p2 + .05*n*p2*randn())
m3 = 200

A2 = zeros(n,n)
A2l = zeros(0,3)
E2 = rand(1:n,m2,2)
for l in 1:size(E2)[1]
	i,j = sort(E2[l,:])
	if i != j && A2[i,j] == 0.
		A2[i,j] = 1.
		A2[j,i] = 1.
global 		A2l = [A2l;[i j 1.];[j i 1.]]
	end
end
@info "Pairwise done."

A3 = zeros(n,n,n)
A3l = zeros(0,4)
A3gt = zeros(0,4)
E3 = rand(1:n,m3,3)
for l in 1:size(E3)[1]
	i,j,k = sort(E3[l,:])
	if i != j && i != k && j != k && A3[i,j,k] == 0.
		A3[i,j,k] = 1.
		A3[i,k,j] = 1.
		A3[j,i,k] = 1.
		A3[j,k,i] = 1.
		A3[k,i,j] = 1.
		A3[k,j,i] = 1.
global		A3l = [A3l;[i j k 1.];[i k j 1.];[j i k 1.];[j k i 1.];[k i j 1.];[k j i 1.]]
global		A3gt = [A3gt;[i j k 1.];[j i k 1.];[k i j 1.]]
	end
end
@info "Triadic done."
# =#

# #= ############### Wheel Hypergraph ##########################
p1 = .3
p2 = .2
p3 = .2
A2l, A3l = gen_rand_hyperwheel_list(n,p1,p2,p3)
A3gt = A3l[[(A3l[i,2] .< A3l[i,3]) for i in 1:size(A3l)[1]],:]
# =#

X = zeros(n,0)
Y = zeros(n,0)

δt = .01
maxiter = 200
for i in 1:100
	x,y = hyper_k(A2l,A3l,zeros(n),.1*rand(n),1.,0.,π/4,π/4,δt,maxiter,1e-6)
global 	X = [X x]
global 	Y = [Y y]
end
@info "Simulations done."


@info "============ STARTING INFERENCE ==========="
ooi = [2,3]
dmax = 2
nkeep = 3000
coeff,ids,e,ee = hyper_inf_filter(X,Y,ooi,dmax,nkeep,.01,.1)

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
tpr3,fpr3,v,I0 = my_ROC_extended(abs.(A3this),A3gt,n)
tpr,fpr,v,I0 = my_ROC_extended(abs.(A2this),A2l,abs.(A3this),A3gt,n)

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




