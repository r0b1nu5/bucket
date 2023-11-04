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
	xxx = hyper_inf(Xs,Ys,ooi,3,1e-1)
	
	
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
	for ij in E2
		i,j = ij
		A2[i,j] = 1.
		A2[j,i] = 1.
	end
	A3 = zeros(n,n,n)
	for ijk in E3
		i,j,k = ijk
		A3[i,j,k] = 1.
		A3[j,k,i] = 1.
		A3[k,i,j] = 1.
		A3[i,k,j] = 1.
		A3[j,i,k] = 1.
		A3[k,j,i] = 1.
	end
	
	
	A2i,A3i = adj_tensors(xxx[1],n,d)
	adji = cat_As(A2i,A3i)
	adj0 = cat_As(A2,A3)
	
	ξ = 1e-10
	roc2 = roc(A2i + ξ*rand(Float64,size(A2i)),A2)
	roc3 = roc(A3i + ξ*rand(Float64,size(A3i)),A3)
	roc23 = roc(adji + ξ*rand(Float64,size(adji)),adj0)

	figure("Hypernetwork of Lorenz oscillators",(14,4))
	subplot(1,3,1)
	PyPlot.plot(roc2.FPR,roc2.TPR,color=cm(l/(n_iter-1)))
	subplot(1,3,2)
	PyPlot.plot(roc3.FPR,roc3.TPR,color=cm(l/(n_iter-1)))
	subplot(1,3,3)
	PyPlot.plot(roc23.FPR,roc23.TPR,color=cm(l/(n_iter-1)))
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



