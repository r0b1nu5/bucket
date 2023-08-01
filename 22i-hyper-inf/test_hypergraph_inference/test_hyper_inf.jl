using Random

include("hyper_inf.jl")
include("hyper_kuramoto.jl")
include("hyper_ktanh.jl")
include("gen_rand_hyperg.jl")
include("tools_hyper.jl")

# ARNI
include("reconstruct.jl")
include("reconstruct_3rd.jl")

# Generating the hypergraph.
n = 7

# Wheel graph with n vertices with randomly added 3-body interactions (p1), and randomly removed 2-edges (p2,p3).
 #=
ntw = "Hyper-wheel"
p1 = .3
p2 = .3 
p3 = .3
# =#

# Wheel graph with randomly removed edges (p2,p3). No 3rd-order edges.
 #=
ntw = "Wheel"
p1 = 0.
p2 = 0.05
p3 = .3
# =# 

# ER random graph with parameter p2. No 3rd-order edges.
 #=
ntw = "ER"
p1 = 0.
p2 = .99
# =#

# ER random graph with parameter p2. Additional 3rd-order edges (proportion p1 of all the possible 3-edges).
# #=
ntw = "Hyper-ER"
p1 = .05
p2 = .4
# =#

if ntw in ["Wheel", "Hyper-wheel"]
	A2,A3 = gen_rand_hyperwheel(n,p1,p2,p3,true)
	A4 = zeros(n,n,n,n)
elseif ntw in ["ER", "Hyper-ER"]
	A2,A3 = gen_hyper_er(n,p1,p2,true)
	A4 = zeros(n,n,n,n)
end

cmapme = get_cmap("RdPu")
cmaparni = get_cmap("GnBu")
cmapdiff = get_cmap("Greys")
cmapme = get_cmap("plasma")
cmaparni = get_cmap("viridis")
#cmapdiff = get_cmap("plasma")
c0 = 0.
c1 = 1.

adj = get_adj_3rd(A2,A3)[1]
# ========================================================================

T = 400
ΔT = 5
NT = 80
# #=
# Testing the efficiency of the inference using the result of the vector field directly.

# Generate the data
# #=
amplitude = .1
X = amplitude*(rand(n,400) .- .5)
Y = f_kuramoto_3rd(X,A2,A3,zeros(n),π/4,π/4)
# =#
 #=
amplitude = .1
ξ0 = 0.0005
X = amplitude*(rand(n,400) .- .5)
Y = f_kuramoto_3rd(X,A2,A3,zeros(n),π/4,π/4) + ξ0*randn(size(Y))
 #=
X = zeros(n,0)
Y = zeros(n,0)
for t in 1:NT
	xxx = hyper_k(A2,A3,zeros(n),amplitude*(rand(n) .- .5),0.,0.,π/4,π/4,.01,ΔT)
	global X = [X xxx[1]]
	global Y = [Y (xxx[2] + ξ0*randn(n,ΔT))]
end
# =#

# =#
 #=
ω = 2*rand(n)
ω .-= mean(ω)
xxx = hyper_k(A2,A3,ω,zeros(n),1.)
X = .2*(rand(n,400) .- .5) + repeat(xxx[1][:,end],1,400)
Y = f_kuramoto_3rd(X,A2,A3,ω)
# =#
 #=
X = amplitude*(rand(n,400) .- .5)
Y = f_ktanh_3rd(X,A2,A3,zeros(n))
# =#
 #=
ω = 2*rand(n)
ω .-= mean(ω)
xxx = hyper_ktanh(A2,A3,ω,zeros(n),1.)
X = .2*(rand(n,400) .- .5) + repeat(xxx[1][:,end],1,400)
Y = f_ktanh_3rd(X,A2,A3,ω)
# =#

test_arni = true
if test_arni
	writedlm("data/test-arni-Xs.csv",X,',')
	writedlm("data/test-arni-Ys.csv",Y,',')
	writedlm("data/test-arni-A2.csv",A2,',')
	writedlm("data/test-arni-A3.csv",A3,',')
end

sen2 = Float64[]
spe2 = Float64[]
sen3 = Float64[]
spe3 = Float64[]
sen4 = Float64[]
spe4 = Float64[]


# Compute the sensitivity and specificity of the inference for various lengths of time series.
iters = 10:5:200
iters = 10:5:80
#iters = 200:5:200
ooi = [2,3]
c = 0
for iter in iters
	global c += 1
	@info "Run $c/$(length(iters))"

	xxx = hyper_inf(X[:,1:iter],Y[:,1:iter],ooi,4,-1e-4)
	A2us = inferred_adj_2nd(xxx[1][2],n)[1]
	A3us = inferred_adj_3rd(xxx[1][3],n)[1]
	adjus = get_adj_3rd(A2us,A3us)[1]

	if test_arni
		adjarni = zeros(n,Int64(n*(n-1)/2))
		for i in 1:n
# bases = ["polynomial", "polynomial_diff", "fourier", "fourier_diff", "power_series", "RBF"]
			w = reconstruct_3rd(X[:,1:iter],Y[:,1:iter],i,adj,1e-6,"power_series")
			adjarni[i,:] = w[1]
		end
		A2arni,A3arni = adj2As(adjarni)
		
		ξ = 1e-10
		
		rocadjus = roc(adjus + ξ*rand(Float64,size(adj)),adj)
		rocA2us = roc(A2us + ξ*rand(n,n),A2)
		rocA3us = roc(A3us + ξ*rand(n,n,n),A3)

		rocadjarni = roc(adjarni + ξ*rand(Float64,size(adj)),adj)
		rocA2arni = roc(A2arni + ξ*rand(n,n),A2)
		rocA3arni = roc(A3arni + ξ*rand(n,n,n),A3)

		figure("ROCs-"*ntw*"-$n")
		subplot(2,3,1)
		PyPlot.plot(rocadjus.FPR,rocadjus.TPR,color=cmapme(c0+c1*iter/maximum(iters)))
		subplot(2,3,2)
		PyPlot.plot(rocA2us.FPR,rocA2us.TPR,color=cmapme(c0+c1*iter/maximum(iters)))
		subplot(2,3,3)
		PyPlot.plot(rocA3us.FPR,rocA3us.TPR,color=cmapme(c0+c1*iter/maximum(iters)))

		subplot(2,3,4)
		PyPlot.plot(rocadjarni.FPR,rocadjarni.TPR,color=cmaparni(c0+c1*iter/maximum(iters)))
		subplot(2,3,5)
		PyPlot.plot(rocA2arni.FPR,rocA2arni.TPR,color=cmaparni(c0+c1*iter/maximum(iters)))
		subplot(2,3,6)
		PyPlot.plot(rocA3arni.FPR,rocA3arni.TPR,color=cmaparni(c0+c1*iter/maximum(iters)))



	end

	yyy = check_inference_bool(A2,A3,A4,xxx[1])
	push!(sen2,yyy[1][1])
	push!(spe2,yyy[1][2])
	push!(sen3,yyy[2][1])
	push!(spe3,yyy[2][2])
	push!(sen4,yyy[3][1])
	push!(spe4,yyy[3][2])
end

figure("ROCs-"*ntw*"-$n")
subplot(2,3,1)
ylabel("TPR")
title("ROC adj, us")
subplot(2,3,2)
title("ROC A2, us")
subplot(2,3,3)
title("ROC A3, us")
subplot(2,3,4)
xlabel("FPR")
ylabel("TPR")
title("ROC adj, ARNI")
subplot(2,3,5)
xlabel("FPR")
title("ROC A2, ARNI")
subplot(2,3,6)
xlabel("FPR")
title("ROC A3, ARNI")

# Plot the sensitivity and specificity as a function of the length of the time series.
#figure("Perfect measurement",(7.5,5))
 #=
PyPlot.plot(iters,sen2,"-o",color="C0",label="2nd-order sen.")
PyPlot.plot(iters,spe2,"--s",color="C0",label="2nd-order spe.")
# =#
 #=
PyPlot.plot(iters,sen3,"-o",color="C1",label="3rd-order sen.")
PyPlot.plot(iters,spe3,"--s",color="C1",label="3rd-order spe.")
# =#
 #=
PyPlot.plot(iters,sen4,"-o",color="C2",label="4th-order sen.")
PyPlot.plot(iters,spe4,"--s",color="C2",label="4th-order spe.")
# =#
 #=
xlabel("Number of measurements")
legend()
# =#

# ========================================================================

 #=
# Testing the efficiency of the inference using the actual time step (extracted from RK4 integration).

# Generate the data (200 time series of length 100 time steps).
X = zeros(n,0)
Y = zeros(n,0)
for i in 1:200
	x,dx = hyper_k(A2,A3,zeros(n),.2*rand(n) .- .1,.1,.0,.001,100)
	global X = [X x]
	global Y = [Y dx]
end

 #=
# Keep one time step for each of the time series.
Xt = X[:,1:100:20000]
Yt = Y[:,1:100:20000]
sens = Float64[]
spes = Float64[]

# Compute the sensitivity and specificity of the inference for various number of time series.
iters = 10:20:200
for iter in iters
	xxx = hyper_inf(Xt[:,1:iter],Yt[:,1:iter],[3,],4)
	yyy = check_inference(A2,A3,xxx[1])
	push!(sens,yyy[2][1])
	push!(spes,yyy[2][2])
end

figure("Actual variation (RK4 measurement)")
PyPlot.plot(iters,sens,"-o",label="sensitivity")
PyPlot.plot(iters,spes,"-o",label="specificity")
xlabel("Number of time steps of length 1")
legend()
# =#

 #=
# Keep two time steps for each of the time series.
idx = sort([1:100:20000;2:100:20000])
Xt = X[:,idx]
Yt = Y[:,idx]
sens = Float64[]
spes = Float64[]

# Compute the sensitivity and specificity of the inference for various number of time series.
iters = 10:20:400
for iter in iters
	xxx = hyper_inf(Xt[:,1:iter],Yt[:,1:iter],[3,],4)
	yyy = check_inference(A2,A3,xxx[1])
	push!(sens,yyy[2][1])
	push!(spes,yyy[2][2])
end

figure("Actual variation (RK4 measurement)")
PyPlot.plot(iters./2,sens,"-o",label="sensitivity")
PyPlot.plot(iters./2,spes,"-o",label="specificity")
xlabel("Number of time steps of length 2")
legend()
# =#

# #=
# Keep five time steps for each of the time series.
idx = sort([1:100:20000;2:100:20000;3:100:20000;4:100:20000;5:100:20000])
Xt = X[:,idx]
Yt = Y[:,idx]
sens = Float64[]
spes = Float64[]

# Compute the sensitivity and specificity of the inference for various number of time series.
iters = 10:50:1000
for iter in iters
	xxx = hyper_inf(Xt[:,1:iter],Yt[:,1:iter],[3,],4)
	yyy = check_inference(A2,A3,xxx[1])
	push!(sens,yyy[2][1])
	push!(spes,yyy[2][2])
end

figure("Actual variation (RK4 measurement)")
PyPlot.plot(iters./5,sens,"-o",label="sensitivity")
PyPlot.plot(iters./5,spes,"-o",label="specificity")
xlabel("Number of time steps of length 5")
legend()
# =#


# Inferred increment: 1) generate a time series with 10x higher resolution and keep 1 out of 10 time steps. 2) infer derivative from time series. 3) run inference.

=#


# ========================================================================

 #=
# Testing the efficiency of the inference using the result of the vector field directly, with fixed time series length and increasing radius of the neighborhood of the fixed point.

sen2 = Float64[]
spe2 = Float64[]
sen3 = Float64[]
spe3 = Float64[]
sen4 = Float64[]
spe4 = Float64[]

# Compute the sensitivity and specificity of the inference for various lengths of time series.
iter = 200
ooi = [3,4]
c = 0

magnitudes = [2e-5,5e-5,1e-4,2e-4,5e-4,1e-3,2e-3,5e-3,1e-2,2e-2,2e-5,.1,.2,.5,.7,.9,1.1,1.3,1.5]

for mag in magnitudes
	global c += 1
	@info "Run $c/$(length(magnitudes))"

	# Generate the data
	X = mag*(rand(n,400) .- .5)
	Y = f_kuramoto_3rd(X,A2,A3,zeros(n))

	xxx = hyper_inf(X[:,1:iter],Y[:,1:iter],ooi,4,1e-6)
	yyy = check_inference_bool(A2,A3,A4,xxx[1])
	push!(sen2,yyy[1][1])
	push!(spe2,yyy[1][2])
	push!(sen3,yyy[2][1])
	push!(spe3,yyy[2][2])
	push!(sen4,yyy[3][1])
	push!(spe4,yyy[3][2])
end

# Plot the sensitivity and specificity as a function of the length of the time series.
figure("Perfect measurement",(7.5,5))
 #=
PyPlot.plot(iters,sen2,"-o",color="C0",label="2nd-order sen.")
PyPlot.plot(iters,spe2,"--s",color="C0",label="2nd-order spe.")
# =#
# #=
PyPlot.semilogx(magnitudes,sen3,"-o",color="C1",label="3rd-order sen.")
PyPlot.plot(magnitudes,spe3,"--s",color="C1",label="3rd-order spe.")
# =#
# #=
PyPlot.plot(magnitudes,sen4,"-o",color="C2",label="4th-order sen.")
PyPlot.plot(magnitudes,spe4,"--s",color="C2",label="4th-order spe.")
# =#
xlabel("Magnitude of the neighborhood")
legend()
# =#


