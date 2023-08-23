using DelimitedFiles, ROC

include("hyper_kuramoto.jl")
include("tools_hyper.jl")
include("hyper_inf.jl")

# ARNI
include("reconstruct_3rd.jl")

BASIS = "power_series"

adj1 = Bool.(readdlm("data/test_ARNI_1.csv",',')) # Only 2nd order interactions
adj2 = Bool.(readdlm("data/test_ARNI_2.csv",',')) # Only 3rd order interactions in triangles
adj3 = Bool.(readdlm("data/test_ARNI_3.csv",',')) # Both 2nd and 3rd order interaction in triangles

n = size(adj1)[1]
m = Int64(n*(n-1)/2)
T = 400
amplitude = .2

A21,A31 = adj2As(adj1)
A22,A32 = adj2As(adj2)
A23,A33 = adj2As(adj3)

X = amplitude*(rand(n,T) .- .5)
Y1 = f_kuramoto_3rd(X,A21,A31,zeros(n),π/4,π/4)
Y2 = f_kuramoto_3rd(X,A22,A32,zeros(n),π/4,π/4)
Y3 = f_kuramoto_3rd(X,A23,A33,zeros(n),π/4,π/4)

# #=
########################## TESTING ARNI WITH 3RD ORDER ####################

adjarni1 = zeros(n,Int64(n*(n-1)/2))
adjarni2 = zeros(n,Int64(n*(n-1)/2))
adjarni3 = zeros(n,Int64(n*(n-1)/2))

th = -1e-4
for i in 1:n
	w1 = reconstruct_3rd(X,Y1,i,adj1,th,BASIS)
	w2 = reconstruct_3rd(X,Y2,i,adj2,th,BASIS)
	w3 = reconstruct_3rd(X,Y3,i,adj3,th,BASIS)

	adjarni1[i,:] = w1[1]
	adjarni2[i,:] = w2[1]
	adjarni3[i,:] = w3[1]
end
A2arni1,A3arni1 = adj2As(adjarni1)
A2arni2,A3arni2 = adj2As(adjarni2)
A2arni3,A3arni3 = adj2As(adjarni3)

ξ = 1e-10
roc11 = roc(adjarni1 .+ ξ*rand(n,m),adj1)
roc12 = roc(adjarni1 .+ ξ*rand(n,m),adj2)
roc13 = roc(adjarni1 .+ ξ*rand(n,m),adj3)
roc21 = roc(adjarni2 .+ ξ*rand(n,m),adj1)
roc22 = roc(adjarni2 .+ ξ*rand(n,m),adj2)
roc23 = roc(adjarni2 .+ ξ*rand(n,m),adj3)
roc31 = roc(adjarni3 .+ ξ*rand(n,m),adj1)
roc32 = roc(adjarni3 .+ ξ*rand(n,m),adj2)
roc33 = roc(adjarni3 .+ ξ*rand(n,m),adj3)

rocA2_11 = roc(A2arni1 + ξ*rand(n,n),A21)
rocA2_12 = roc(A2arni1 + ξ*rand(n,n),A22)
rocA2_13 = roc(A2arni1 + ξ*rand(n,n),A23)
rocA2_21 = roc(A2arni2 + ξ*rand(n,n),A21)
rocA2_22 = roc(A2arni2 + ξ*rand(n,n),A22)
rocA2_23 = roc(A2arni2 + ξ*rand(n,n),A23)
rocA2_31 = roc(A2arni3 + ξ*rand(n,n),A21)
rocA2_32 = roc(A2arni3 + ξ*rand(n,n),A22)
rocA2_33 = roc(A2arni3 + ξ*rand(n,n),A23)

rocA3_11 = roc(A3arni1 + ξ*rand(n,n,n),A31)
rocA3_12 = roc(A3arni1 + ξ*rand(n,n,n),A32)
rocA3_13 = roc(A3arni1 + ξ*rand(n,n,n),A33)
rocA3_21 = roc(A3arni2 + ξ*rand(n,n,n),A31)
rocA3_22 = roc(A3arni2 + ξ*rand(n,n,n),A32)
rocA3_23 = roc(A3arni2 + ξ*rand(n,n,n),A33)
rocA3_31 = roc(A3arni3 + ξ*rand(n,n,n),A31)
rocA3_32 = roc(A3arni3 + ξ*rand(n,n,n),A32)
rocA3_33 = roc(A3arni3 + ξ*rand(n,n,n),A33)

figure("ARNI w/ 3rd order")
subplot(3,3,1)
PyPlot.plot(roc11.FPR,roc11.TPR,label="AUC = $(AUC(roc11))")
PyPlot.plot(rocA2_11.FPR,rocA2_11.TPR,label="AUC-A2 = $(AUC(rocA2_11))")
PyPlot.plot(rocA3_11.FPR,rocA3_11.TPR,label="AUC-A3 = $(AUC(rocA3_11))")
legend()
title("ARNI 1 vs. adj 1")
subplot(3,3,2)
PyPlot.plot(roc12.FPR,roc12.TPR,label="AUC = $(AUC(roc12))")
PyPlot.plot(rocA2_12.FPR,rocA2_12.TPR,label="AUC-A2 = $(AUC(rocA2_12))")
PyPlot.plot(rocA3_12.FPR,rocA3_12.TPR,label="AUC-A3 = $(AUC(rocA3_12))")
legend()
title("ARNI 1 vs. adj 2")
subplot(3,3,3)
PyPlot.plot(roc13.FPR,roc13.TPR,label="AUC = $(AUC(roc13))")
PyPlot.plot(rocA2_13.FPR,rocA2_13.TPR,label="AUC-A2 = $(AUC(rocA2_13))")
PyPlot.plot(rocA3_13.FPR,rocA3_13.TPR,label="AUC-A3 = $(AUC(rocA3_13))")
legend()
title("ARNI 1 vs. adj 3")
subplot(3,3,4)
PyPlot.plot(roc21.FPR,roc21.TPR,label="AUC = $(AUC(roc21))")
PyPlot.plot(rocA2_21.FPR,rocA2_21.TPR,label="AUC-A2 = $(AUC(rocA2_21))")
PyPlot.plot(rocA3_21.FPR,rocA3_21.TPR,label="AUC-A3 = $(AUC(rocA3_21))")
legend()
title("ARNI 2 vs. adj 1")
subplot(3,3,5)
PyPlot.plot(roc22.FPR,roc22.TPR,label="AUC = $(AUC(roc22))")
PyPlot.plot(rocA2_22.FPR,rocA2_22.TPR,label="AUC-A2 = $(AUC(rocA2_22))")
PyPlot.plot(rocA3_22.FPR,rocA3_22.TPR,label="AUC-A3 = $(AUC(rocA3_22))")
legend()
title("ARNI 2 vs. adj 2")
subplot(3,3,6)
PyPlot.plot(roc23.FPR,roc23.TPR,label="AUC = $(AUC(roc23))")
PyPlot.plot(rocA2_23.FPR,rocA2_23.TPR,label="AUC-A2 = $(AUC(rocA2_23))")
PyPlot.plot(rocA3_23.FPR,rocA3_23.TPR,label="AUC-A3 = $(AUC(rocA3_23))")
legend()
title("ARNI 2 vs. adj 3")
subplot(3,3,7)
PyPlot.plot(roc31.FPR,roc31.TPR,label="AUC = $(AUC(roc31))")
PyPlot.plot(rocA2_31.FPR,rocA2_31.TPR,label="AUC-A2 = $(AUC(rocA2_31))")
PyPlot.plot(rocA3_31.FPR,rocA3_31.TPR,label="AUC-A3 = $(AUC(rocA3_31))")
legend()
title("ARNI 3 vs. adj 1")
subplot(3,3,8)
PyPlot.plot(roc32.FPR,roc32.TPR,label="AUC = $(AUC(roc32))")
PyPlot.plot(rocA2_32.FPR,rocA2_32.TPR,label="AUC-A2 = $(AUC(rocA2_32))")
PyPlot.plot(rocA3_32.FPR,rocA3_32.TPR,label="AUC-A3 = $(AUC(rocA3_32))")
legend()
title("ARNI 3 vs. adj 2")
subplot(3,3,9)
PyPlot.plot(roc33.FPR,roc33.TPR,label="AUC = $(AUC(roc33))")
PyPlot.plot(rocA2_33.FPR,rocA2_33.TPR,label="AUC-A2 = $(AUC(rocA2_33))")
PyPlot.plot(rocA3_33.FPR,rocA3_33.TPR,label="AUC-A3 = $(AUC(rocA3_33))")
legend()
title("ARNI 3 vs. adj 3")


# =#

##################### OUR METHOD #############################
# #=

adjus1 = zeros(n,Int64(n*(n-1)/2))
adjus2 = zeros(n,Int64(n*(n-1)/2))
adjus3 = zeros(n,Int64(n*(n-1)/2))

th = -1e-4
ooi = [2,3]
w1 = hyper_inf(X,Y1,ooi,3,-1e-4)
w2 = hyper_inf(X,Y2,ooi,3,-1e-4)
w3 = hyper_inf(X,Y3,ooi,3,-1e-4)

A2us1 = inferred_adj_2nd(w1[1][2],n)[1]
A3us1 = inferred_adj_3rd(w1[1][3],n)[1]
A2us2 = inferred_adj_2nd(w2[1][2],n)[1]
A3us2 = inferred_adj_3rd(w2[1][3],n)[1]
A2us3 = inferred_adj_2nd(w3[1][2],n)[1]
A3us3 = inferred_adj_3rd(w3[1][3],n)[1]

adjus1 = get_adj_3rd(A2us1,A3us1)[1]
adjus2 = get_adj_3rd(A2us2,A3us2)[1]
adjus3 = get_adj_3rd(A2us3,A3us3)[1]


ξ = 1e-10
roc11 = roc(adjus1 .+ ξ*rand(n,m),adj1)
roc12 = roc(adjus1 .+ ξ*rand(n,m),adj2)
roc13 = roc(adjus1 .+ ξ*rand(n,m),adj3)
roc21 = roc(adjus2 .+ ξ*rand(n,m),adj1)
roc22 = roc(adjus2 .+ ξ*rand(n,m),adj2)
roc23 = roc(adjus2 .+ ξ*rand(n,m),adj3)
roc31 = roc(adjus3 .+ ξ*rand(n,m),adj1)
roc32 = roc(adjus3 .+ ξ*rand(n,m),adj2)
roc33 = roc(adjus3 .+ ξ*rand(n,m),adj3)

rocA2_11 = roc(A2us1 + ξ*rand(n,n),A21)
rocA2_12 = roc(A2us1 + ξ*rand(n,n),A22)
rocA2_13 = roc(A2us1 + ξ*rand(n,n),A23)
rocA2_21 = roc(A2us2 + ξ*rand(n,n),A21)
rocA2_22 = roc(A2us2 + ξ*rand(n,n),A22)
rocA2_23 = roc(A2us2 + ξ*rand(n,n),A23)
rocA2_31 = roc(A2us3 + ξ*rand(n,n),A21)
rocA2_32 = roc(A2us3 + ξ*rand(n,n),A22)
rocA2_33 = roc(A2us3 + ξ*rand(n,n),A23)

rocA3_11 = roc(A3us1 + ξ*rand(n,n,n),A31)
rocA3_12 = roc(A3us1 + ξ*rand(n,n,n),A32)
rocA3_13 = roc(A3us1 + ξ*rand(n,n,n),A33)
rocA3_21 = roc(A3us2 + ξ*rand(n,n,n),A31)
rocA3_22 = roc(A3us2 + ξ*rand(n,n,n),A32)
rocA3_23 = roc(A3us2 + ξ*rand(n,n,n),A33)
rocA3_31 = roc(A3us3 + ξ*rand(n,n,n),A31)
rocA3_32 = roc(A3us3 + ξ*rand(n,n,n),A32)
rocA3_33 = roc(A3us3 + ξ*rand(n,n,n),A33)

figure("US w/ 3rd order")
subplot(3,3,1)
PyPlot.plot(roc11.FPR,roc11.TPR,label="AUC = $(AUC(roc11))")
PyPlot.plot(rocA2_11.FPR,rocA2_11.TPR,label="AUC-A2 = $(AUC(rocA2_11))")
PyPlot.plot(rocA3_11.FPR,rocA3_11.TPR,label="AUC-A3 = $(AUC(rocA3_11))")
legend()
title("us 1 vs. adj 1")
subplot(3,3,2)
PyPlot.plot(roc12.FPR,roc12.TPR,label="AUC = $(AUC(roc12))")
PyPlot.plot(rocA2_12.FPR,rocA2_12.TPR,label="AUC-A2 = $(AUC(rocA2_12))")
PyPlot.plot(rocA3_12.FPR,rocA3_12.TPR,label="AUC-A3 = $(AUC(rocA3_12))")
legend()
title("us 1 vs. adj 2")
subplot(3,3,3)
PyPlot.plot(roc13.FPR,roc13.TPR,label="AUC = $(AUC(roc13))")
PyPlot.plot(rocA2_13.FPR,rocA2_13.TPR,label="AUC-A2 = $(AUC(rocA2_13))")
PyPlot.plot(rocA3_13.FPR,rocA3_13.TPR,label="AUC-A3 = $(AUC(rocA3_13))")
legend()
title("us 1 vs. adj 3")
subplot(3,3,4)
PyPlot.plot(roc21.FPR,roc21.TPR,label="AUC = $(AUC(roc21))")
PyPlot.plot(rocA2_21.FPR,rocA2_21.TPR,label="AUC-A2 = $(AUC(rocA2_21))")
PyPlot.plot(rocA3_21.FPR,rocA3_21.TPR,label="AUC-A3 = $(AUC(rocA3_21))")
legend()
title("us 2 vs. adj 1")
subplot(3,3,5)
PyPlot.plot(roc22.FPR,roc22.TPR,label="AUC = $(AUC(roc22))")
PyPlot.plot(rocA2_22.FPR,rocA2_22.TPR,label="AUC-A2 = $(AUC(rocA2_22))")
PyPlot.plot(rocA3_22.FPR,rocA3_22.TPR,label="AUC-A3 = $(AUC(rocA3_22))")
legend()
title("us 2 vs. adj 2")
subplot(3,3,6)
PyPlot.plot(roc23.FPR,roc23.TPR,label="AUC = $(AUC(roc23))")
PyPlot.plot(rocA2_23.FPR,rocA2_23.TPR,label="AUC-A2 = $(AUC(rocA2_23))")
PyPlot.plot(rocA3_23.FPR,rocA3_23.TPR,label="AUC-A3 = $(AUC(rocA3_23))")
legend()
title("us 2 vs. adj 3")
subplot(3,3,7)
PyPlot.plot(roc31.FPR,roc31.TPR,label="AUC = $(AUC(roc31))")
PyPlot.plot(rocA2_31.FPR,rocA2_31.TPR,label="AUC-A2 = $(AUC(rocA2_31))")
PyPlot.plot(rocA3_31.FPR,rocA3_31.TPR,label="AUC-A3 = $(AUC(rocA3_31))")
legend()
title("us 3 vs. adj 1")
subplot(3,3,8)
PyPlot.plot(roc32.FPR,roc32.TPR,label="AUC = $(AUC(roc32))")
PyPlot.plot(rocA2_32.FPR,rocA2_32.TPR,label="AUC-A2 = $(AUC(rocA2_32))")
PyPlot.plot(rocA3_32.FPR,rocA3_32.TPR,label="AUC-A3 = $(AUC(rocA3_32))")
legend()
title("us 3 vs. adj 2")
subplot(3,3,9)
PyPlot.plot(roc33.FPR,roc33.TPR,label="AUC = $(AUC(roc33))")
PyPlot.plot(rocA2_33.FPR,rocA2_33.TPR,label="AUC-A2 = $(AUC(rocA2_33))")
PyPlot.plot(rocA3_33.FPR,rocA3_33.TPR,label="AUC-A3 = $(AUC(rocA3_33))")
legend()
title("us 3 vs. adj 3")

# =#


