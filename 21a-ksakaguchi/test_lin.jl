include("iterations.jl")

ntw = "ntw10"
#ntw = "ntw9"

L = readdlm("ntw_data/"*ntw*"_L.csv",',')
include("ntw_data/"*ntw*"_cycles.jl")
ω = vec(readdlm("ntw_data/"*ntw*"_om1.csv",','))
#ω = rand(10)
#ω .-= mean(ω)
x0 = 2π*rand(10)

b,w = L2B(L) # Undirected incidence matrix
B,Bout,Bin = L2B_bidir(L) # Directed incidence matrices
Ld = pinv(L)

n,m = size(B)
n,m2 = size(b)

a1 = 40*rand(m2)
a2 = rand(m2)
A = diagm(0 => [a1;a2])
La = Bout*A*B'
xf = Ld*ω

Δf = b'*xf # Real angle differences
Ff = A*[Δf;-Δf] # Real flows

#u = [0,1,0] # Tentative winding vector
T = 1000 # Number of iterations
#T = 20
λ = .1 # Scaling parameter for the fixed point iteration Tu

#Lmin = 1.0*ones(m2) # Scaling for the cutset projection. Lower bounds on the coupling derivatives.
Lmin = 1.0*ones(m)

ρ1 = .1
δ = .1

#Δ01 = (2*rand(m2) .- 1)*γ
Δ01 = 4π*rand(m2) .- 2π
Δ1 = Array{Float64,2}(undef,m2,0)
Δ1 = [Δ1 b'*Ld*b*Δ01]
for t in 1:T
	global Δ1
	Δ = Δ1[:,end] - δ*b'*Ld*(Bout*A*[Δ1[:,end];-Δ1[:,end]] - ω)
	Δ1 = [Δ1 Δ]
end
f1 = A*[Δ1;-Δ1]

#Δ02 = (2*rand(m2) .- 1)*γ
Δ02 = 4π*rand(m2) .- 2π
Δ4 = Array{Float64,2}(undef,m2,0)
Δ4 = [Δ4 b'*Ld*b*Δ02]
for t in 1:T
	global Δ4
	Δ = Δ4[:,end] - δ*b'*Ld*(Bout*A*[Δ4[:,end];-Δ4[:,end]] - ω)
	Δ4 = [Δ4 Δ]
end
f4 = A*[Δ4;-Δ4]

##=
#t1,t2,t3 = check_monotonicity(f1,f4,P)
#@info "$t1, $t2, $t3"

df1 = [norm(f1[:,t]-Ff) for t in 1:T]
df4 = [norm(f4[:,t]-Ff) for t in 1:T]
dΔ1 = [norm(dcc(Δ1[:,t]-Δf)) for t in 1:T]
dΔ4 = [norm(dcc(Δ4[:,t]-Δf)) for t in 1:T]
# #=
df14 = [norm(f1[:,t]-f4[:,t]) for t in 1:T]
δΔ14 = [norm(Δ1[:,t]-Δ4[:,t]) for t in 1:T]

df14i = [norm(f1[:,t]-f4[:,t],Inf) for t in 1:T]
δΔ14i = [norm(Δ1[:,t]-Δ4[:,t],Inf) for t in 1:T]
# =#


figure(ntw,(12.,8.))
subplot(2,3,1)
PyPlot.semilogy((1:T),df1,label="||f - f*||","o",color="C0")
PyPlot.plot((1:T),df4,label="||f - f*||","x",color="C1")

legend()
xlabel("iteration")
ylabel("error")

subplot(2,3,2)
PyPlot.semilogy((1:T),df14,label="||f1 - f2||_2","-o")
legend()
xlabel("iteration")

subplot(2,3,3)
PyPlot.semilogy((1:T),df14i,label="||f1 - f2||_∞","-o")
legend()
xlabel("iteration")

subplot(2,3,4)
PyPlot.semilogy((1:T),dΔ1,label="||Δ - Δ*||","o",color="C2")
PyPlot.plot((1:T),dΔ4,label="||Δ - Δ*||","x",color="C3")

legend()
xlabel("iteration")
ylabel("error")

subplot(2,3,5)
PyPlot.semilogy((1:T),δΔ14,label="||Δ1 - Δ2||_2","-o",color="C3")
legend()
xlabel("iteration")

subplot(2,3,6)
PyPlot.semilogy((1:T),δΔ14i,label="||Δ1 - Δ2||_∞","-o",color="C3")
legend()
xlabel("iteration")
# =#
