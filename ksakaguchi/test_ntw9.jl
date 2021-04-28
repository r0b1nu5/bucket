include("ksakaguchi.jl")
include("iterations.jl")

L = readdlm("ntw_data/ntw9_L.csv",',')
include("ntw_data/ntw9_cycles.jl")
ω = rand(9)
ω .-= mean(ω)
θ0 = 2π*rand(9)

b,w = L2B(L) # Undirected incidence matrix
B,Bout,Bin = L2B_bidir(L) # Directed incidence matrices

n,m = size(B)
n,m2 = size(b)

α = .1 # Phase frustration
γ = π/2 - α # Bound on the angle difference

xxx = ksakaguchi(L,ω,θ0,α)
θf = xxx[1][:,end] # Real solution
Δf = dcc(b'*θf) # Real angle differences
Ff = h([Δf;-Δf]) # Real flows
uf = winding(θf,Σ) # Real winding vector

#u = [0,1,0] # Tentative winding vector
u = uf
T = 200 # Number of iterations
λ = 1. # Scaling parameter for the fixed point iteration Tu

Lmin = 1.0*ones(m2) # Scaling for the cutset projection. Lower bounds on the coupling derivatives.

f01 = 2*γ*rand(m) .- γ
f1,f2,Δ1,Δ2 = iterations(f01,Bout,b,ω,u,Lmin,γ,λ,T)

f02 = 2*γ*rand(m) .- γ
f3,f4,Δ3,Δ4 = iterations(f02,Bout,b,ω,u,Lmin,γ,λ,T)

df01 = norm(f01 - Ff)
df1 = [norm(f1[:,t]-Ff) for t in 1:T]
df2 = [norm(f2[:,t]-Ff) for t in 1:T]
dΔ1 = [norm(Δ1[:,t] - [Δf;-Δf]) for t in 1:T]
dΔ2 = [norm(Δ2[:,t] - Δf) for t in 1:T]

df02 = norm(f02 - Ff)
df3 = [norm(f3[:,t]-Ff) for t in 1:T]
df4 = [norm(f4[:,t]-Ff) for t in 1:T]
dΔ3 = [norm(Δ3[:,t] - [Δf;-Δf]) for t in 1:T]
dΔ4 = [norm(Δ4[:,t] - Δf) for t in 1:T]

df00 = norm(f01 - f02)
df13 = [norm(f1[:,t]-f3[:,t]) for t in 1:T]
df24 = [norm(f2[:,t]-f4[:,t]) for t in 1:T]
dΔ13 = [norm(Δ1[:,t]-Δ3[:,t]) for t in 1:T]
dΔ24 = [norm(Δ2[:,t]-Δ4[:,t]) for t in 1:T]

figure()
subplot(1,3,1)
PyPlot.semilogy((0:T),[df01;df1],label="||f - f*||")
PyPlot.plot((0:T),[df01;df2],label="||f' - f*||")
PyPlot.plot((1:T),dΔ1,label="||Δ - Δ*||")
PyPlot.plot((1:T),dΔ2,label="||Δ' - Δ*||")
legend()
xlabel("iteration")
ylabel("error")

subplot(1,3,2)
PyPlot.semilogy((0:T),[df02;df3],label="||f - f*||")
PyPlot.plot((0:T),[df02;df4],label="||f' - f*||")
PyPlot.plot((1:T),dΔ3,label="||Δ - Δ*||")
PyPlot.plot((1:T),dΔ4,label="||Δ' - Δ*||")
legend()
xlabel("iteration")

subplot(1,3,3)
PyPlot.semilogy((0:T),[df00;df13],label="||f1 - f2||")
PyPlot.plot((0:T),[df00;df24],label="||f1' - f2'||")
PyPlot.plot((1:T),dΔ13,label="||Δ1 - Δ2||")
PyPlot.plot((1:T),dΔ24,label="||Δ1' - Δ2'||")
legend()
xlabel("iteration")


