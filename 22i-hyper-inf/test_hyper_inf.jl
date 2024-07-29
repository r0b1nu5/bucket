using Random

include("hyper_inf.jl")
include("hyper_kuramoto.jl")
include("hyper_ktanh.jl")
include("gen_rand_hyperg.jl")
include("tools_hyper.jl")

include("arni-reconstruct.jl")
include("arni-reconstruct-3rd.jl")

# Generating the hypergraph.
n = 10
 #=
ntw = "Hyper-wheel"
p1 = .3
p2 = .3 
p3 = .3
# =#
 #=
ntw = "Wheel"
p1 = 0.
p2 = 0.05
p3 = .3
# =# 
 #=
ntw = "ER"
p1 = 0.
p2 = .99
# =#
# #=
ntw = "Hyper-ER"
p1 = .05
p2 = .4
# =#
 #= 
ntw = "Hyper-ER"
p1 = .5
p2 = .8
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
cmaparni = get_cmap("plasma")
#cmapdiff = get_cmap("plasma")
c0 = 0.
c1 = 1.

#adj = get_adj_3rd(A2,A3)[1]
adj = cat_As(A2,A3)
# ========================================================================

# #=
# Testing the efficiency of the inference using the result of the vector field directly.

# Generate the data
 #=
amplitude = .1
X = amplitude*(rand(n,400) .- .5)
Y = f_kuramoto_3rd(X,A2,A3,zeros(n),π/4,π/4) + .01*randn(size(X))
# =#
# #= ########## PERFECT MEASUREMENTS ###############
amplitude = 2.
ξ0 = 0.0005
X = amplitude*(rand(n,400) .- .5)
Y = f_kuramoto_3rd(X,A2,A3,zeros(n),π/4,π/4) + ξ0*randn(size(X))
# =#
 #= ########### TRUNCATE TIME SERIES IN THE δ-BOX #############
T = 10000
δ0 = .1
δ1 = 2.
h = .001
amplitude = 2π
X = zeros(n,0)
Yh = zeros(n,0)
Y = zeros(n,0)
count = 0
while size(X)[2] < T
	global count += 1
	xxx = hyper_k(2*A2,A3,zeros(n),amplitude*(rand(n) .- .5),0.,0.,π/4,π/4,h,Int64(50/h))
	xx = mod.(xxx[1] .+ π,2π) .- π
	xx = xx - repeat(xx[:,end],1,size(xx)[2])
	t = [(δ0 < norm(xx[:,i]) < δ1) for i in 1:(size(xx)[2]-1)]
	idx = vec(1:size(xx)[2]-1)[t]
	yyy = (xxx[1][:,2:end] - xxx[1][:,1:end-1])./h
	global X = [X xx[:,idx]]
	global Yh = [Yh yyy[:,idx]]
	global Y = [Y xxx[2][:,idx]]
end
@info "Number of time series used: $count"
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

test_arni = false
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
#iters = 10:5:200
#iters = 10:15:150
#iters = 5000:1000:T
#iters = 10:10:200
iters = 200:5:200
ooi = [2,3]
dmax = 2
c = 0
for iter in iters
	global c += 1
	@info "Run $c/$(length(iters))"

	xxx = hyper_inf(X[:,1:iter],Y[:,1:iter],ooi,dmax,1e-1)
#	A2us = inferred_adj_2nd(xxx[1][2],n)[1]
#	A2us = inferred_adj_2nd(xxx[1][2],n)[2]
	A2us = xxx[1][2]
#	A3us = inferred_adj_3rd(xxx[1][3],n)[1]
#	A3us = inferred_adj_3rd(xxx[1][3],n)[2]
	A3us = xxx[1][3]
	adjus = get_adj_3rd(A2us,A3us)[1]
	adju = cat_As(A2us,A3us)
@info "============= WE ARE DONE ================"

	if test_arni
		adjarni = zeros(n,Int64(n*(n-1)/2))
		for i in 1:n
			w = reconstruct_3rd(X[:,1:iter],Y[:,1:iter],i,adj,1e-6,"power_series")
#			w = reconstruct_3rd(X[:,1:2*iter],Y[:,1:2*iter],i,adj,1e-6,"power_series")
			adjarni[i,:] = w[1]
		end
		A2arni,A3arni = adj2As(adjarni)
		adja = cat_As(A2arni,A3arni)

		ξ = 1e-10
		
#		rocadjus = roc(adjus + ξ*rand(Float64,size(adj)),adj)
		rocadjus = roc(abs.(adju) + ξ*rand(Float64,size(adju)),adj)
		rocA2us = roc(abs.(A2us) + ξ*rand(n,n),A2)
		rocA3us = roc(abs.(A3us) + ξ*rand(n,n,n),A3)

#		rocadjarni = roc(adjarni + ξ*rand(Float64,size(adj)),adj)
		rocadjarni = roc(adja + ξ*rand(Float64,size(adj)),adj)
		rocA2arni = roc(A2arni + ξ*rand(n,n),A2)
		rocA3arni = roc(A3arni + ξ*rand(n,n,n),A3)

		figure("ROCs-"*ntw*"-$n",(15,10))
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
	else
		ξ = 1e-10
		
#		rocadjus = roc(adjus + ξ*rand(Float64,size(adj)),adj)
		rocadjus = roc(abs.(adju) + ξ*rand(Float64,size(adju)),adj)
		rocA2us = roc(abs.(A2us) + ξ*rand(n,n),A2)
		rocA3us = roc(abs.(A3us) + ξ*rand(n,n,n),A3)

		figure("ROCs-"*ntw*"-$n",(15,4))
		subplot(1,3,1)
		PyPlot.plot(rocadjus.FPR,rocadjus.TPR,color=cmapme(c0+c1*iter/maximum(iters)))
		subplot(1,3,2)
		PyPlot.plot(rocA2us.FPR,rocA2us.TPR,color=cmapme(c0+c1*iter/maximum(iters)))
		subplot(1,3,3)
		PyPlot.plot(rocA3us.FPR,rocA3us.TPR,color=cmapme(c0+c1*iter/maximum(iters)))

	end

	yyy = check_inference_bool(A2,A3,A4,xxx[1])
	push!(sen2,yyy[1][1])
	push!(spe2,yyy[1][2])
	push!(sen3,yyy[2][1])
	push!(spe3,yyy[2][2])
	push!(sen4,yyy[3][1])
	push!(spe4,yyy[3][2])
end

if test_arni
	figure("ROCs-"*ntw*"-$n",(15,10))
	subplot(2,3,1)
	ylabel("TPR")
	title("ROC adj, THIS")
	subplot(2,3,2)
	title("ROC A2, THIS")
	subplot(2,3,3)
	title("ROC A3, THIS")
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
else
	figure("ROCs-"*ntw*"-$n",(15,10))
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
end

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

	xxx = hyper_inf(X[:,1:iter],Y[:,1:iter],ooi,4,1e-1)
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


