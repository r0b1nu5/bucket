include("hyper_inf.jl")
include("hyper_kuramoto.jl")
include("gen_rand_hyperg.jl")

n = 7
p1 = .5
p2 = .2 
p3 = .2

A2,A3 = gen_rand_hyperwheel(n,p1,p2,p3)

 #=
# Testing the efficiency of the inference using the result of the vector field directly.

X = .2*rand(n,100) .- .1
Y = f_kuramoto_3rd(X,A2,A3,zeros(n))

sens = Float64[]
spes = Float64[]
iters = 10:5:100
for iter in iters
	xxx = hyper_inf(X[:,1:iter],Y[:,1:iter],[3,],4)
	yyy = check_inference(A2,A3,xxx[1])
	push!(sens,yyy[2][1])
	push!(spes,yyy[2][2])
end

figure("Perfect measurement")
PyPlot.plot(iters,sens,"-o",label="sensitivity")
PyPlot.plot(iters,spes,"-o",label="specificity")
xlabel("Number of measurements")
legend()
# =#

# #=
# Testing the efficiency of the inference using the actual time step (extracted from RK4 integration).

X = zeros(n,0)
Y = zeros(n,0)

for i in 1:200
	x,dx = hyper_k(A2,A3,zeros(n),.2*rand(n) .- .1,.1,.0,.001,100)
	global X = [X x]
	global Y = [Y dx]
end

 #=
Xt = X[:,1:100:20000]
Yt = Y[:,1:100:20000]
sens = Float64[]
spes = Float64[]
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
idx = sort([1:100:20000;2:100:20000])
Xt = X[:,idx]
Yt = Y[:,idx]
sens = Float64[]
spes = Float64[]
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
idx = sort([1:100:20000;2:100:20000;3:100:20000;4:100:20000;5:100:20000])
Xt = X[:,idx]
Yt = Y[:,idx]
sens = Float64[]
spes = Float64[]
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






