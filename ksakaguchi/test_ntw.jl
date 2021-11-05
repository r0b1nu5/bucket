include("ksakaguchi.jl")
include("iterations.jl")

ntw = "ntw10"
#ntw = "ntw9"

L = readdlm("ntw_data/"*ntw*"_L.csv",',')
include("ntw_data/"*ntw*"_cycles.jl")
ω = vec(readdlm("ntw_data/"*ntw*"_om1.csv",','))
#ω = rand(10)
#ω .-= mean(ω)
θ0 = vec(readdlm("ntw_data/"*ntw*"_th1.csv",','))
#θ0 = 2π*rand(10)

b,w = L2B(L) # Undirected incidence matrix
B,Bout,Bin = L2B_bidir(L) # Directed incidence matrices

n,m = size(B)
n,m2 = size(b)

α = .1 # Phase frustration
γ = π/2 - α # Bound on the angle difference
hγ1 = h(-γ,α)
hγ2 = h(γ,α)

xxx = ksakaguchi(L,ω,θ0,α,false,false,.01,1e-8)
θf = xxx[1][:,end] # Real solution
uf = winding(θf,Σ) # Real winding vector

Δf = dcc(b'*θf) # Real angle differences
Ff = h([Δf;-Δf]) # Real flows

#u = [0,1,0] # Tentative winding vector
u = uf
T = 1000 # Number of iterations
#T = 20
λ = .1 # Scaling parameter for the fixed point iteration Tu

#Lmin = 1.0*ones(m2) # Scaling for the cutset projection. Lower bounds on the coupling derivatives.
Lmin = 1.0*ones(m)
P = dir_cycle_proj(B,Lmin)

 #=
f01 = (hγ2 - hγ1)*rand(m) .+ hγ1
f01 = (-1. + sin(α))*ones(m) + .1*rand(m)
f1,f2,Δ1,Δ2 = iterations1(f01,Bout,b,ω,u,Lmin,γ,λ,T)

f02 = (hγ2 - hγ1)*rand(m) .+ hγ1
f02 = (-1. + sin(α))*ones(m) + .1*rand(m)
f3,f4,Δ3,Δ4 = iterations1(f02,Bout,b,ω,u,Lmin,γ,λ,T)
# =#

 #=
Δ01 = 2*γ*rand(Int(m/2)) .- γ
Δ01 = -2*γ .+ .1*rand(m2)
f1,f2,Δ1,Δ2 = iterations2(Δ01,Bout,b,ω,u,Lmin,γ,λ,T)

Δ02 = 2*γ*rand(Int(m/2)) .- γ
Δ02 = -2*γ .+ .1*rand(m2)
f3,f4,Δ3,Δ4 = iterations2(Δ02,Bout,b,ω,u,Lmin,γ,λ,T)

df1 = [norm(f1[:,t]-Ff) for t in 1:T]
df2 = [norm(f2[:,t]-Ff) for t in 1:T]
df3 = [norm(f3[:,t]-Ff) for t in 1:T]
df4 = [norm(f4[:,t]-Ff) for t in 1:T]

dΔ1 = [norm(Δ1[:,t] - [Δf;-Δf]) for t in 1:T]
dΔ2 = [norm(Δ2[:,t] - Δf) for t in 1:T]
dΔ3 = [norm(Δ3[:,t] - [Δf;-Δf]) for t in 1:T]
dΔ4 = [norm(Δ4[:,t] - Δf) for t in 1:T]

df13 = [norm(f1[:,t]-f3[:,t]) for t in 1:T]
df24 = [norm(f2[:,t]-f4[:,t]) for t in 1:T]
dΔ13 = [norm(Δ1[:,t]-Δ3[:,t]) for t in 1:T]
dΔ24 = [norm(Δ2[:,t]-Δ4[:,t]) for t in 1:T]

df13i = [norm(f1[:,t]-f3[:,t],Inf) for t in 1:T]
df24i = [norm(f2[:,t]-f4[:,t],Inf) for t in 1:T]
dΔ13i = [norm(Δ1[:,t]-Δ3[:,t],Inf) for t in 1:T]
dΔ24i = [norm(Δ2[:,t]-Δ4[:,t],Inf) for t in 1:T]

figure()
subplot(1,4,1)
PyPlot.semilogy((1:T),df1,label="||f - f*||","o")
PyPlot.plot((1:T),df2,label="||f' - f*||","o")
PyPlot.plot((1:T),dΔ1,label="||Δ - Δ*||","o")
PyPlot.plot((1:T),dΔ2,label="||Δ' - Δ*||","o")
legend()
xlabel("iteration")
ylabel("error")

subplot(1,4,2)
PyPlot.semilogy((1:T),df3,label="||f - f*||","o")
PyPlot.plot((1:T),df4,label="||f' - f*||","o")
PyPlot.plot((1:T),dΔ3,label="||Δ - Δ*||","o")
PyPlot.plot((1:T),dΔ4,label="||Δ' - Δ*||","o")
legend()
xlabel("iteration")

subplot(1,4,3)
PyPlot.semilogy((1:T),df13,label="||f1 - f2||_2","-o")
PyPlot.plot((1:T),df24,label="||f1' - f2'||_2","-o")
PyPlot.plot((1:T),dΔ13,label="||Δ1 - Δ2||_2","-o")
PyPlot.plot((1:T),dΔ24,label="||Δ1' - Δ2'||_2","-o")
legend()
xlabel("iteration")

subplot(1,4,4)
PyPlot.semilogy((1:T),df13i,label="||f1 - f2||_∞","-o")
PyPlot.plot((1:T),df24i,label="||f1' - f2'||_∞","-o")
PyPlot.plot((1:T),dΔ13i,label="||Δ1 - Δ2||_∞","-o")
PyPlot.plot((1:T),dΔ24i,label="||Δ1' - Δ2'||_∞","-o")
legend()
xlabel("iteration")

# =#


# #=
#f01 = (hγ2 - hγ1)*rand(m) .+ hγ1
#f01 = hγ1 .+ .01*rand(m)
#f01 = [hγ1 .+ .01*rand(m2);hγ2 .- .01*rand(m2)]
#f1,f2,f3 = iterations3(f01,Bout,B,ω,u,Lmin,γ,λ,T)
#f1 = iterations4(f01,θf,Bout,B,ω,u,Lmin,γ,λ,T,true)

ρ1 = .1
ρ2 = .1

#Δ01 = (2*rand(m2) .- 1)*γ
Δ01 = 4π*rand(m2) .- 2π
Δ1 = iterations5(Δ01,θf,Bout,B,C,ω,α,u,ρ1,T,false)
f1 = h([Δ1;-Δ1])

#f02 = (hγ2 - hγ1)*rand(m) .+ hγ1
#f02 = hγ1 .+ .01*rand(m)
#f02 = [hγ1 .+ .01*rand(m2);hγ2 .- .01*rand(m2)]
#f4,f5,f6 = iterations3(f02,Bout,B,ω,u,Lmin,γ,λ,T)
#f4 = iterations4(f02,θf,Bout,B,ω,u,Lmin,γ,λ,T,false)

#Δ02 = (2*rand(m2) .- 1)*γ
Δ02 = 4π*rand(m2) .- 2π
Δ4 = iterations5(Δ02,θf,Bout,B,C,ω,α,u,ρ1,T,false)
f4 = h([Δ4;-Δ4])

##=
#t1,t2,t3 = check_monotonicity(f1,f4,P)
#@info "$t1, $t2, $t3"

df1 = [norm(f1[:,t]-Ff) for t in 1:T]
df4 = [norm(f4[:,t]-Ff) for t in 1:T]
dΔ1 = [norm(dcc(Δ1[:,t]-Δf)) for t in 1:T]
dΔ4 = [norm(dcc(Δ4[:,t]-Δf)) for t in 1:T]
# #=
df14 = [norm(f1[:,t]-f4[:,t]) for t in 1:T]
dΔ14 = [norm(dcc(Δ1[:,t]-Δ4[:,t])) for t in 1:T]
δΔ14 = [norm(Δ1[:,t]-Δ4[:,t]) for t in 1:T]

df14i = [norm(f1[:,t]-f4[:,t],Inf) for t in 1:T]
dΔ14i = [norm(dcc(Δ1[:,t]-Δ4[:,t]),Inf) for t in 1:T]
δΔ14i = [norm(Δ1[:,t]-Δ4[:,t],Inf) for t in 1:T]
# =#

 #=
P = dir_cycle_proj(B,ones(m))
R = pinv(P)*P

df14 = [norm(R*(f1[:,t]-f4[:,t])) for t in 1:T]
df25 = [norm(R*(f2[:,t]-f5[:,t])) for t in 1:T]
df36 = [norm(R*(f3[:,t]-f6[:,t])) for t in 1:T]

df14i = [norm(R*(f1[:,t]-f4[:,t]),Inf) for t in 1:T]
df25i = [norm(R*(f2[:,t]-f5[:,t]),Inf) for t in 1:T]
df36i = [norm(R*(f3[:,t]-f6[:,t]),Inf) for t in 1:T]
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
PyPlot.semilogy((1:T),dΔ14,label="||dcc(Δ1 - Δ2)||_2","-o",color="C2")
PyPlot.semilogy((1:T),δΔ14,label="||Δ1 - Δ2||_2","-o",color="C3")
legend()
xlabel("iteration")

subplot(2,3,6)
PyPlot.semilogy((1:T),dΔ14i,label="||dcc(Δ1 - Δ2)||_∞","-o",color="C2")
PyPlot.semilogy((1:T),δΔ14i,label="||Δ1 - Δ2||_∞","-o",color="C3")
legend()
xlabel("iteration")
# =#
