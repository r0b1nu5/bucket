include("hyper_inf.jl")
include("hyper_kuramoto.jl")
include("hyper_ktanh.jl")
include("gen_rand_hyperg.jl")
include("tools_hyper.jl")

include("../../ARNI/reconstruct.jl")

# Generating the hypergraph.
n = 7
 #=
p1 = .3
p2 = .3 
p3 = .3
# =#
# #=
p1 = 0.
p2 = 0.05
p3 = .3
# =# 
amplitude = .2

test_arni = true
if test_arni
	figure("ROCs")
	subplot(2,3,1)
	title("us")
	xlabel("FPR")
	ylabel("TPR")
	subplot(2,3,2)
	title("ARNI (polynomial)")
	xlabel("FPR")
	ylabel("TPR")
	subplot(2,3,3)
	title("ARNI (polynomial diff)")
	xlabel("FPR")
	ylabel("TPR")
	subplot(2,3,4)
	title("ARNI (fourier)")
	xlabel("FPR")
	ylabel("TPR")
	subplot(2,3,5)
	title("ARNI (fourier diff)")
	xlabel("FPR")
	ylabel("TPR")
	subplot(2,3,6)
	title("ARNI (poly) - us")
	xlabel("FPR")
	ylabel("TPR (ARNI-us)")
end

cmapme = get_cmap("RdPu")
cmaparni = get_cmap("GnBu")
cmapdiff = get_cmap("Greys")
cmapme = get_cmap("plasma")
cmaparni = get_cmap("viridis")
#cmapdiff = get_cmap("plasma")
c0 = 0.
c1 = 1.

# #=
A2,A3 = gen_rand_hyperwheel(n,p1,p2,p3,true)
A4 = zeros(n,n,n,n)
# =#
 #=
A2 = diagm(0 => n*ones(n)) - ones(n,n)
A3 = zeros(n,n,n)
A4 = zeros(n,n,n,n)
# =#

adj = get_adj_3rd(A2,A3)
# ========================================================================

# #=
# Testing the efficiency of the inference using the result of the vector field directly.

# Generate the data
 #=
X = amplitude*(rand(n,400) .- .5)
Y = f_kuramoto_3rd(X,A2,A3,zeros(n))
# =#
# #=
X = amplitude*(rand(n,400) .- .5)
Y = f_kuramoto_3rd(X,A2,A3,zeros(n),π/4,π/4)
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

if test_arni
	writedlm("data/test-arni-Xs.csv",X,',')
	writedlm("data/test-arni-Ys.csv",Y,',')
	writedlm("data/test-arni-adj.csv",A2,',')
end

sen2 = Float64[]
spe2 = Float64[]
sen3 = Float64[]
spe3 = Float64[]
sen4 = Float64[]
spe4 = Float64[]

fprarni1 = zeros(n^2+1,0)
tprarni1 = zeros(n^2+1,0)
aucarni1 = Float64[]
fprarni2 = zeros(n^2+1,0)
tprarni2 = zeros(n^2+1,0)
aucarni2 = Float64[]
fprarni3 = zeros(n^2+1,0)
tprarni3 = zeros(n^2+1,0)
aucarni3 = Float64[]
fprarni4 = zeros(n^2+1,0)
tprarni4 = zeros(n^2+1,0)
aucarni4 = Float64[]
fprme = zeros(n^2+1,0)
tprme = zeros(n^2+1,0)
aucme = Float64[]

# Compute the sensitivity and specificity of the inference for various lengths of time series.
iters = 10:5:200
iters = 10:5:80
ooi = [2,]
c = 0
for iter in iters
	global c += 1
	@info "Run $c/$(length(iters))"

	xxx = hyper_inf(X[:,1:iter],Y[:,1:iter],ooi,4,-1e-4)
	if test_arni
		 #=
		adjme = zeros(n,n)
		for k in keys(xxx[1][2])
			adjme[k[2][1],k[2][2]] = xxx[1][2][k]
		end
		# =#
		# #=
	    	adjme = inferred_adj_2nd(xxx[1][2],n)[1]
		# =#
		
		adjarni1 = zeros(n,n)
		adjarni2 = zeros(n,n)
		adjarni3 = zeros(n,n)
		adjarni4 = zeros(n,n)
		for i in 1:n
# bases = ["polynomial", "polynomial_diff", "fourier", "fourier_diff", "power_series", "RBF"]
			w = reconstruct(X[:,1:iter],Y[:,1:iter],i,A2,"polynomial")
			#w = reconstruct(X[:,1:iter],Y[:,1:iter],i,A2,"power_series")
			adjarni1[i,:] = w[1]
			w = reconstruct(X[:,1:iter],Y[:,1:iter],i,A2,"polynomial_diff")
			adjarni2[i,:] = w[1]
			w = reconstruct(X[:,1:iter],Y[:,1:iter],i,A2,"fourier")
			adjarni3[i,:] = w[1]
			w = reconstruct(X[:,1:iter],Y[:,1:iter],i,A2,"fourier_diff")
			adjarni4[i,:] = w[1]
		end
		
		ξ = 1e-10
		rocme = roc(adjme .+ ξ*rand(n,n),A2) # Adding rand(n,n) makes sure that FPR and TPR have the same dim as fprme and tprme.
		x0 = rocme.FPR
		y0 = rocme.TPR
		global fprme = [fprme rocme.FPR]
		global tprme = [tprme rocme.TPR]
		push!(aucme,AUC(rocme))
		
		rocarni = roc(adjarni1 .+ ξ*rand(n,n),A2)
		x1 = rocarni.FPR
		y1 = rocarni.TPR
		global fprarni1 = [fprarni1 rocarni.FPR]
		global tprarni1 = [tprarni1 rocarni.TPR]
		push!(aucarni1,AUC(rocarni))
		
		rocarni = roc(adjarni2 .+ ξ*rand(n,n),A2)
		x2 = rocarni.FPR
		y2 = rocarni.TPR
		global fprarni2 = [fprarni2 rocarni.FPR]
		global tprarni2 = [tprarni2 rocarni.TPR]
		push!(aucarni2,AUC(rocarni))
		
		rocarni = roc(adjarni3 .+ ξ*rand(n,n),A2)
		x3 = rocarni.FPR
		y3 = rocarni.TPR
		global fprarni3 = [fprarni3 rocarni.FPR]
		global tprarni3 = [tprarni3 rocarni.TPR]
		push!(aucarni3,AUC(rocarni))
		
		rocarni = roc(adjarni4 .+ ξ*rand(n,n),A2)
		x4 = rocarni.FPR
		y4 = rocarni.TPR
		global fprarni4 = [fprarni4 rocarni.FPR]
		global tprarni4 = [tprarni4 rocarni.TPR]
		push!(aucarni4,AUC(rocarni))
		

		Z = sortslices([x1 y1 -ones(length(y1));x0 -ones(length(y0)) y0],dims=1)
		x = Z[:,1]
		z1 = Z[:,2]
		temp = Int64[]
		for i in 1:length(z1)
			if z1[i] < 0
				push!(temp,i)
			else
				z1[temp] .= z1[i]
				temp = Int64[]
			end
		end
		z1[temp] .= 1.
		z0 = Z[:,3]
		temp = Int64[]
		for i in 1:length(z0)
			if z0[i] < 0
				push!(temp,i)
			else
				z0[temp] .= z0[i]
				temp = Int64[]
			end
		end
		z0[temp] .= 1.

		figure("ROCs")
		subplot(2,3,1)
		PyPlot.plot(rocme.FPR,rocme.TPR,color=cmapme(c0 + c1*iter/maximum(iters)))
		subplot(2,3,2)
		PyPlot.plot(fprarni1[:,end],tprarni1[:,end],color=cmaparni(c0 + c1*iter/maximum(iters)))
		subplot(2,3,3)
		PyPlot.plot(fprarni2[:,end],tprarni2[:,end],color=cmaparni(c0 + c1*iter/maximum(iters)))
		subplot(2,3,4)
		PyPlot.plot(fprarni3[:,end],tprarni3[:,end],color=cmaparni(c0 + c1*iter/maximum(iters)))
		subplot(2,3,5)
		PyPlot.plot(fprarni4[:,end],tprarni4[:,end],color=cmaparni(c0 + c1*iter/maximum(iters)))
		subplot(2,3,6)
		PyPlot.plot(x,z0-z1,color=cmapdiff(1-.9*(c0 + c1*iter/maximum(iters))))
	end




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
# #=
PyPlot.plot(iters,sen2,"-o",color="C0",label="2nd-order sen.")
PyPlot.plot(iters,spe2,"--s",color="C0",label="2nd-order spe.")
# =#
# #=
PyPlot.plot(iters,sen3,"-o",color="C1",label="3rd-order sen.")
PyPlot.plot(iters,spe3,"--s",color="C1",label="3rd-order spe.")
# =#
# #=
PyPlot.plot(iters,sen4,"-o",color="C2",label="4th-order sen.")
PyPlot.plot(iters,spe4,"--s",color="C2",label="4th-order spe.")
# =#
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


