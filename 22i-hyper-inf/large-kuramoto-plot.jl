using DelimitedFiles, PyPlot, Distributed, Dates


#n = 300; run = "008"
#n = 300; run = "011"
#n = 500; run = "009"
#n = 300; run = "013"
n = 300; run = "014"
nks = [2000,5000,10000]

cmap = get_cmap("viridis")

figure("ROK",(13.5,4))
subplot(1,3,1)
PyPlot.plot([0,1],[0,1],"--k",lw=.5)
subplot(1,3,2)
PyPlot.plot([0,1],[0,1],"--k",lw=.5)
subplot(1,3,3)
PyPlot.plot([0,1],[0,1],"--k",lw=.5)

c = length(nks)+1
for nkeep in nks[end:-1:1]
	tpr2 = readdlm("data/rok-"*run*"-keep$nkeep-tpr2.csv",',')
	fpr2 = readdlm("data/rok-"*run*"-keep$nkeep-fpr2.csv",',')
	tpr3 = readdlm("data/rok-"*run*"-keep$nkeep-tpr3.csv",',')
	fpr3 = readdlm("data/rok-"*run*"-keep$nkeep-fpr3.csv",',')
	tpr = readdlm("data/rok-"*run*"-keep$nkeep-tpr.csv",',')
	fpr = readdlm("data/rok-"*run*"-keep$nkeep-fpr.csv",',')

	global c -= 1
	col = cmap(c/(length(nks)+1))

	figure("ROK")
	subplot(1,3,1)
	PyPlot.plot(fpr,tpr,color=col,lw=2,label="nk = $nkeep")
	PyPlot.plot([fpr[end],1],[tpr[end],1],"--",color=col,lw=2)
	subplot(1,3,2)
	PyPlot.plot(fpr2,tpr2,color=col)
	PyPlot.plot([fpr2[end],1],[tpr2[end],1],"--",color=col)
	subplot(1,3,3)
	PyPlot.plot(fpr3,tpr3,color=col)
	PyPlot.plot([fpr3[end],1],[tpr3[end],1],"--",color=col)
end

subplot(1,3,1)
xlabel("FPR")
ylabel("TPR")
legend()
title("ROC: all edges")
subplot(1,3,2)
xlabel("FPR")
title("ROC: pairwise")
subplot(1,3,3)
xlabel("FPR")
title("ROC: triadic")

PyPlot.savefig("temp/rok-"*run*"-$(now()).pdf",format="pdf")


