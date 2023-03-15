#include("hyper_inf.jl")
include("hyper_inf_new.jl")
include("hyper_kuramoto.jl")
include("gen_rand_hyperg.jl")

# Generating the hypergraph.
n = 7
# #=
p1 = .3
p2 = .3 
p3 = .3
# =#
 #=
p1 = 0.
p2 = 0.
p3 = .3
# =# 

A2,A3 = gen_rand_hyperwheel(n,p1,p2,p3,true)
A4 = zeros(n,n,n,n)

# ========================================================================

# #=
# Testing the efficiency of the inference using the result of the vector field directly.

# Generate the data
X = π*(rand(n,400) .- .5)  
Y = f_kuramoto_3rd(X,A2,A3,zeros(n))

sen2 = Float64[]
spe2 = Float64[]
sen3 = Float64[]
spe3 = Float64[]
sen4 = Float64[]
spe4 = Float64[]

# Compute the sensitivity and specificity of the inference for various lengths of time series.
iters = 10:5:200
ooi = [2,3]
c = 0
for iter in iters
	global c += 1
	@info "Run $c/$(length(iters))"

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


