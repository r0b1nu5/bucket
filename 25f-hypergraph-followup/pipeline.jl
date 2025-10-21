using PyPlot

include("gen-hypergraph.jl")
include("hyper-kuramoto.jl")
include("tools.jl")


n = 10
p = .2
ω0 = .0
ξ0 = .5

h = .01

A2 = zeros(0,3)
A3,B,B2,E = rand_3_digraph(n,p)
ω1 = ω0*rand(n); ω1 .-= mean(ω1)
ϕ2 = .2
ϕ3 = .2

# 0. Reach the sync state
Θ1,dΘ1,iter1 = hyper_k(A2,A3,ω1,zeros(n),ϕ2,ϕ3,h,10000,1e-5)
θstar = Θ1[:,end]

# 1. Run the system without control
ω2 = ω1 .- mean(dΘ1[:,end])
Θ2,dΘ2,iter2 = hyper_k_gaussian_noise(A2,A3,ω2,Θ1[:,end],ξ0,ϕ2,ϕ3,h,10000,-1.)

for i in 1:n
	PyPlot.plot(h*(1:length(Θ2[i,:])),Θ2[i,:],color="C$(i-1)")
end

# 2. Run THIS
b = ones(n)
d = .8*ones(n)

idx = get_loose_ends(A2,A3)
b = zeros(n)
d = zeros(n)
if length(idx) == 0
	@info "System is connected"
	b[1] = 1.
	d[1] = .8
else
	@info "System is not connected"
	b[idx] = ones(length(idx)) 
	d[idx] = .8*ones(length(idx))
end

# 3. Run the damped system
Θ3,dΘ3,iter3 = hyper_k_drooped_gaussian_noise(A2,A3,ω2,Θ2[:,end],10*b,θstar,ξ0,ϕ2,ϕ3,h,10000,-1.)
#Θ3,dΘ3,iter3 = hyper_k_damped_gaussian_noise(A2,A3,ω2,Θ2[:,end],d,ξ0,ϕ2,ϕ3,h,1000,-1.)

for i in 1:n
	PyPlot.plot(h*(length(Θ2[i,:]) .+ (1:length(Θ3[i,:]))),Θ3[i,:],color="C$(i-1)")
end





