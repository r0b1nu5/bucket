using DelimitedFiles

include("cyclic_iterations_repo.jl")

L = readdlm("ntw_data/ntw11_L.csv",',')
include("ntw_data/ntw11_cycles.jl")
ω = vec(readdlm("ntw_data/ntw11_om1.csv",','))

B,w,Bt = L2B(L)

n,m = size(B)

ϕ = .4*rand() .+ 1
a = .4*rand() .+ .8
h1 = (x -> a*(sin(x - ϕ) + sin(ϕ)))
γ1 = (-π/2 + ϕ,π/2 - ϕ)

ϕs = .4*rand(2*m) .+ .1
as = .4*rand(2*m) .+ .8
h2,hi2,γ2 = load_ksakaguchi(as,ϕs)

Δ1 = (γ1[2] - γ1[1])*rand(m) .+ γ1[1]
u = [0,1,0,0]
δ = .01
s1 = 1.
xxx1 = iterations(Δ1,B,C,u,ω,h1,γ1,δ,s1,1000,1e-6,false)
@info "Uniform coupling: OK!"

Δ2 = [(γ2[i][2] - γ2[i][1])*rand() + γ2[i][1] for i in 1:m]
u = [0,1,0,0]
δ = .01
s2 = ones(2*m)
xxx2 = iterations(Δ2,B,C,u,ω,h2,γ2,δ,s2,1000,1e-6,false)
@info "Heterogeneous coupling: OK!"



