using PyPlot, DelimitedFiles

include("hyper_inf.jl")

n = 100
iters = [3000,4000,5000,6000,7000]

cmapme = get_cmap("plasma")

for iter in iters
	A2 = readdlm("data/kuramoto-n$n-iter$iter-A2.csv",',')
	A2l = zeros(0,3)
	for i in 1:n
		for j in 1:n
			if abs(A2[i,j]) > 1e-10
				A2l = [A2l;[i j A2[i,j]]]
			end
		end
	end
	A2us = readdlm("data/kuramoto-n$n-iter$iter-A2this.csv",',')
	a3 = readdlm("data/kuramoto-n$n-iter$iter-A3.csv",',')
	A3 = zeros(n,n,n)
	for i in 1:n^3
		A3[i] = a3[i]
	end
	A3l = zeros(0,4)
	for i in 1:n
		for j in 1:n-1
			for k in j+1:n
				if i != j && i != k && abs(A3[i,j,k]) > 1e-10
					A3l = [A3l;[i j k A3[i,j,k]]]
				end
			end
		end
	end
	A3us = readdlm("data/kuramoto-n$n-iter$iter-A3this.csv",',')

	ξ = 1e-10

	tpr,fpr = my_ROC(abs.(A2us),A2l,abs.(A3us),A3l,n)
	tpr2,fpr2 = my_ROC(abs.(A2us),A2l,n)
	tpr3,fpr3 = my_ROC(abs.(A3us),A3l,n)

	figure("ROCs-$n",(15,4))
	subplot(1,3,1)
	PyPlot.plot(fpr,tpr,color=cmapme((iter-minimum(iters))/max(1,(maximum(iters)-minimum(iters)))))
	subplot(1,3,2)
	PyPlot.plot(fpr2,tpr2,color=cmapme((iter-minimum(iters))/max(1,(maximum(iters)-minimum(iters)))))
	subplot(1,3,3)
	PyPlot.plot(fpr3,tpr3,color=cmapme((iter-minimum(iters))/max(1,(maximum(iters)-minimum(iters)))))
end


figure("ROCs-$n",(15,10))
subplot(1,3,1)
xlabel("FPR")
ylabel("TPR")
title("ROC adj, THIS")
subplot(1,3,2)
xlabel("FPR")
title("ROC A2, THIS")
subplot(1,3,3)
title("ROC A3, THIS")
xlabel("FPR")
ylabel("TPR")


