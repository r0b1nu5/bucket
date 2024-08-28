using PyPlot, DelimitedFiles, ROC, Dates

include("hyper_inf.jl")

cm = get_cmap("cividis")
n_iter = 1

t = .01
d = 3

λ = .01 # SINDy's sparsity parameter
α = .9 # correlation threshold

t_ratio = Float64[]

 #=
for l in 0:9
XXs = Matrix(readdlm("data/coupled_lorenz_solution-$l.txt")')
x = readdlm("data/lorenz-edges-$l.txt",'}')
T0 = 0
dT = 150
T = 300
n_iter = 10
# =#
 #=
for l in 0:0
XXs = Matrix(readdlm("data/coupled_lorenz_solution_10.txt")')
x = readdlm("data/lorenz_edges_10.txt",'}')
T0 = 100
dT = 700
T = 1000
# =#
# #=
for l in 0:0
N = 10
M = 500
m = 100
tt = 3
T0 = 50
dT = 150
T = 300
XXs = Matrix(readdlm("data/coupled_lorenz_solution_$(N)_$(M)_$(tt).txt")')
x = readdlm("data/lorenz_edges_$(N)_$(M)_$(tt).txt",'}')
# =#

 #=
Xs = XXs[:,1:end-1]
Ys = (XXs[:,2:end]-XXs[:,1:end-1])/t
# =#
# #=
Xs = zeros(size(XXs)[1],0)
Ys = zeros(size(XXs)[1],0)
for k in 1:m
	Xs = [Xs XXs[:,T0 .+ (1:dT) .+ (k-1)*T]]
	Ys = [Ys (XXs[:,T0 .+ (2:dT+1) .+ (k-1)*T] - XXs[:,T0 .+ (1:dT) .+ (k-1)*T])/t]
end
# =#

N,T = size(Ys)
n = Int64(N/d)
	
ooi = [2,3]

@info "Inference starts..."
t0 = time()
xxx = hyper_inf(Xs,Ys,ooi,3,λ)
t1 = time()
t2 = time()
yyy = hyper_inf_filter(Xs,Ys,ooi,3,α,λ)
t3 = time()
push!(t_ratio,(t3-t2)/(t1-t0))

@info "THIS without preprocessing: time = $(t1-t0)'', relative error = $(xxx[4])"
@info "THIS with preprocessing: time = $(t3-t2)'', relative error = $(yyy[4])"

a2us = xxx[1][2]
a = zeros(n,n)
for l in 1:size(a2us)[1]
	i,j = ceil.(Int64,a2us[l,1:2]./d)
	a[i,j] = max(a[i,j],abs(a2us[l,3]))
end
b2us = yyy[1][2]
b = zeros(n,n)
for l in 1:size(b2us)[1]
	i,j = ceil.(Int64,b2us[l,1:2]./d)
	b[i,j] = max(b[i,j],abs(b2us[l,3]))
end

A2us = zeros(0,3)
B2us = zeros(0,3)
for i in 1:n
	for j in 1:n
		if i != j
			A2us = [A2us;[i j a[i,j]]]
			B2us = [B2us;[i j b[i,j]]]
		end
	end
end

a3us = xxx[1][3]
a = zeros(n,n,n)
for l in 1:size(a3us)[1]
	i,j,k = ceil.(Int64,a3us[l,1:3]./d)
	a[i,j,k] = max(a[i,j,k],abs(a3us[l,4]))
end
b3us = yyy[1][3]
b = zeros(n,n,n)
for l in 1:size(b3us)[1]
	i,j,k = ceil.(Int64,b3us[l,1:3]./d)
	b[i,j,k] = max(b[i,j,k],abs(b3us[l,4]))
end

A3us = zeros(0,4)
B3us = zeros(0,4)
for i in 1:n
	for j in 1:n-1
		for k in j+1:n
			if i != j && i != k
				A3us = [A3us;[i j k a[i,j,k]]]
				B3us = [B3us;[i j k b[i,j,k]]]
			end
		end
	end
end
	
#A2,AA2 = inferred_adj_2nd(xxx[1][2],n)
#A3,AA3 = inferred_adj_3rd(xxx[1][3],n)

	
	
E2 = Vector{Int64}[]
for s in x[1,:]
	global e = Int64[]
	for c in s
		if !(c in [',','[',']','{',' '])
			push!(e,parse(Int64,c)+1)
		end
	end
	if length(e) > 0
		push!(E2,e)
	end
end
E3 = Vector{Int64}[]
for s in x[2,:]
	global e = Int64[]
	for c in s
		if !(c in [',','[',']','{',' '])
			push!(e,parse(Int64,c)+1)
		end
	end
	if length(e) > 0
		push!(E3,e)
	end
end

A2 = zeros(n,n)
A2l = zeros(0,3)
for ij in E2
	i,j = sort(ij)
	A2[i,j] = 1.
	A2[j,i] = 1.
	A2l = [A2l;[i j 1.];[j i 1.]]
end
A3 = zeros(n,n,n)
A3l = zeros(0,4)
for ijk in E3
	i,j,k = sort(ijk)
	A3[i,j,k] = 1.
	A3[j,k,i] = 1.
	A3[k,i,j] = 1.
	A3[i,k,j] = 1.
	A3[j,i,k] = 1.
	A3[k,j,i] = 1.
	A3l = [A3l;[i j k 1.];[j i k 1.];[k i j 1.]]
end

tpr,fpr = my_ROC(abs.(A2us),A2l,abs.(A3us),A3l,n)
tpr2,fpr2 = my_ROC(abs.(A2us),A2l,n)
tpr3,fpr3 = my_ROC(abs.(A3us),A3l,n)

tpr_,fpr_ = my_ROC(abs.(B2us),A2l,abs.(B3us),A3l,n)
tpr2_,fpr2_ = my_ROC(abs.(B2us),A2l,n)
tpr3_,fpr3_ = my_ROC(abs.(B3us),A3l,n)

figure("ROCs-Lorenz-10",(15,4))
subplot(2,3,1)
PyPlot.plot(fpr,tpr,color=cm((l+1)/n_iter))
subplot(2,3,2)
PyPlot.plot(fpr2,tpr2,color=cm((l+1)/n_iter))
subplot(2,3,3)
PyPlot.plot(fpr3,tpr3,color=cm((l+1)/n_iter))
	
subplot(2,3,4)
PyPlot.plot(fpr_,tpr_,color=cm((l+1)/n_iter))
subplot(2,3,5)
PyPlot.plot(fpr2_,tpr2_,color=cm((l+1)/n_iter))
subplot(2,3,6)
PyPlot.plot(fpr3_,tpr3_,color=cm((l+1)/n_iter))

end

subplot(2,3,1)
title("ROC: 2-edges")
ylabel("TPR")
subplot(2,3,2)
title("ROC: 3-edges")
subplot(2,3,3)
title("ROC: all edges")
subplot(2,3,4)
xlabel("FPR")
subplot(2,3,5)
xlabel("FPR")
subplot(2,3,6)
xlabel("FPR")




