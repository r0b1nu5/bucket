using DelimitedFiles, PyPlot, Distributed, Dates


#n = 300; run = "008"
#n = 300; run = "011"
n = 500; run = "009"

tpr2 = readdlm("data/rok-"*run*"-tpr2.csv",',')
fpr2 = readdlm("data/rok-"*run*"-fpr2.csv",',')
tpr3 = readdlm("data/rok-"*run*"-tpr3.csv",',')
fpr3 = readdlm("data/rok-"*run*"-fpr3.csv",',')
tpr = readdlm("data/rok-"*run*"-tpr.csv",',')
fpr = readdlm("data/rok-"*run*"-fpr.csv",',')


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



