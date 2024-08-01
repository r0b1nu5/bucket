using PyPlot, DelimitedFiles

include("hyper_inf.jl")
include("hyper_kuramoto.jl")

n = 100
iter = 7000

A2 = readdlm("data/kuramoto-Simplicial-ER-n$n-iter$iter-A2.csv",',')
A2l = zeros(0,3)
for i in 1:n
	for j in 1:n
		if abs(A2[i,j]) > 1e-10
			global A2l = [A2l;[i j A2[i,j]]]
		end
	end
end
A2us = readdlm("data/kuramoto-Simplicial-ER-n$n-iter$iter-A2this.csv",',')

a3 = readdlm("data/kuramoto-Simplicial-ER-n$n-iter$iter-A3.csv",',')
A3 = zeros(n,n,n)
for i in 1:n^3
	A3[i] = a3[i]
end
A3l = zeros(0,4)
for i in 1:n
	for j in 1:n-1
		for k in j+1:n
			if i != j && i != k && abs(A3[i,j,k]) > 1e-10
				global A3l = [A3l;[i j k A3[i,j,k]]]
			end
		end
	end
end
A3us = readdlm("data/kuramoto-Simplicial-ER-n$n-iter$iter-A3this.csv",',')

tpr,fpr,v,I0 = my_ROC_extended(abs.(A2us),A2l,abs.(A3us),A3l,n)


# #= ################# 1st strategy ####################################################
m,id = findmax(tpr-fpr)
vmax = v[id]

figure()
PyPlot.plot(v,tpr,label="TPR")
PyPlot.plot(v,fpr,label="FPR")
PyPlot.plot(v,tpr-fpr,label="TPR-FPR")
PyPlot.plot([vmax,vmax],[0,1],"--k")
xlabel("threshold")
legend()

N1 = sum(abs.(A2us[:,3]) .> vmax) + sum(abs.(A3us[:,4]) .> vmax)
N2 = sum(abs.(A2us[:,3]) .> 1e-10) + sum(abs.(A3us[:,4]) .> 1e-10)
N3 = n*(n-1) + n*binomial(n-1,2)

@info "$N1 hyperedges kept out of $N2 ($(100*round(N1/N2,digits=2))%)."

# =#
 #= ################# 2nd strategy ####################################################
compute_new_ts = false
amplitude = 1.
ξ0 = 0.0005
@info "Computing time series..."
if compute_new_ts
	X = amplitude*(rand(n,iter) .- .5)
	Y = f_kuramoto_3rd(X,A2,A3,zeros(n),π/4,π/4) + ξ0*randn(size(X))
end
@info "Time series computed."

Yh = zeros(size(Y))
no = sum(Y.^2)
re = [1.,]
v = v[abs.(v) .< Inf]
m_max = 400000
dm = 10000
for m in dm:dm:m_max
	@info "Hyperedge no. $m"
	for mm in 1:dm
		e = I0[m-dm+mm]
		Yh[e[1],:] += v[m-dm+mm]*vec(prod(X[e[2:end],:],dims=1))
	end
	push!(re,sum((Y-Yh).^2)/no)
end

figure()
PyPlot.plot(0:dm:m_max,re)
xlabel("number of hyperedges")
ylabel("relative error")

# =#





