using PyPlot, DelimitedFiles, ROC

include("hyper_inf.jl")

cm = get_cmap("cividis")

t = .01
d = 3

n_iter = 10
for l in 0:n_iter-1
	XXs = Matrix(readdlm("data/coupled_lorenz_solution-$l.txt")')
#	XXs = Matrix(readdlm("data/coupled_lorenz_solution-1.txt")')

	Xs = zeros(size(XXs)[1],0)
	Ys = zeros(size(XXs)[1],0)
	for k in 1:10
		Xs = [Xs XXs[:,(1:150) .+ (k-1)*300]]
		Ys = [Ys (XXs[:,(2:151) .+ (k-1)*300] - XXs[:,(1:150) .+ (k-1)*300])/t]
	end

	N,T = size(Ys)
	n = Int64(N/d)
	
	ooi = [2,3]

	@info "Inference starts..."
#	xxx = hyper_inf(Xs,Ys,ooi,3,1e-1)
	xxx = hyper_inf_filter(Xs,Ys,ooi,3,.8,1e-1)

	a2us = xxx[1][2]
	a = zeros(n,n)
	for l in 1:size(a2us)[1]
		i,j = ceil.(Int64,a2us[l,1:2]./d)
		a[i,j] = max(a[i,j],abs(a2us[l,3]))
	end
	A2us = zeros(0,3)
	for i in 1:n
		for j in 1:n
			if i != j
				A2us = [A2us;[i j a[i,j]]]
			end
		end
	end
	a3us = xxx[1][3]
	a = zeros(n,n,n)
	for l in 1:size(a3us)[1]
		i,j,k = ceil.(Int64,a3us[l,1:3]./d)
		a[i,j,k] = max(a[i,j,k],abs(a3us[l,4]))
	end
	A3us = zeros(0,4)
	for i in 1:n
		for j in 1:n-1
			for k in j+1:n
				if i != j && i != k
					A3us = [A3us;[i j k a[i,j,k]]]
				end
			end
		end
	end
	
	#A2,AA2 = inferred_adj_2nd(xxx[1][2],n)
	#A3,AA3 = inferred_adj_3rd(xxx[1][3],n)
	
	
	
	x = readdlm("data/lorenz-edges-$l.txt",'}')
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
	
	ξ = 1e-10

	tpr,fpr = my_ROC(abs.(A2us),A2l,abs.(A3us),A3l,n)
	tpr2,fpr2 = my_ROC(abs.(A2us),A2l,n)
	tpr3,fpr3 = my_ROC(abs.(A3us),A3l,n)

	figure("ROCs-Lorenz-",(15,4))
	subplot(1,3,1)
	PyPlot.plot(fpr,tpr,color=cm((l+1)/n_iter))
	subplot(1,3,2)
	PyPlot.plot(fpr2,tpr2,color=cm((l+1)/n_iter))
	subplot(1,3,3)
	PyPlot.plot(fpr3,tpr3,color=cm((l+1)/n_iter))
	
end
	
subplot(1,3,1)
title("ROC: 2-edges")
xlabel("FPR")
ylabel("TPR")
subplot(1,3,2)
title("ROC: 3-edges")
xlabel("FPR")
subplot(1,3,3)
title("ROC: all edges")
xlabel("FPR")



