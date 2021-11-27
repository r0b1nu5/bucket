using DelimitedFiles

include("acyclic_algorithm_repo.jl")

h = ((x,a) -> (sin(x - a) + sin(a)))
hi = ((f,a) -> (asin(f - sin(a)) + a))

ω = .6*rand(8)
ω .- mean(ω)

L1 = readdlm("ntw_data/ntw8_acyc_L.csv",',')
B1,w1 = L2B(L1) 

IJV = readdlm("ntw_data/ntw8_acyc_Lsp.csv",',')
L2 = sparse(Int64.(IJV[:,1]),Int64.(IJV[:,2]),IJV[:,3])
B2,w2 = L2B(L2)

n,m = size(B1)

α = .3
γ1 = π/2 - α
h1 = ((x -> h(x,α)),(f -> hi(f,α)))
γ1 = (-γ1,γ1)

αs = .4*rand(14) .+ .1
γ2 = [(-π/2 + αs[i],π/2 - αs[i]) for i in 1:2*m]
h2 = Vector{Tuple{Function,Function}}()
for i in 1:m
	push!(h2,((x -> (sin(x - αs[i]) + sin(αs[i]))),(f -> (asin(f - sin(αs[i])) + αs[i]))))
	push!(γ2,(-π/2 + max(αs[i],αs[i+m]),π/2 - max(αs[i],αs[i+m])))
end
for i in m+1:2*m
	push!(h2,((x -> (sin(x - αs[i]) + sin(αs[i]))),(f -> (asin(f - sin(αs[i])) + αs[i]))))
	push!(γ2,(-π/2 + max(αs[i],αs[i-m]),π/2 - max(αs[i],αs[i-m])))
end


xxx1 = run_acyclic_algorithm(B1,ω,h1,γ1)
@info "Dense matrix and uniform coupling: OK!"

xxx2 = run_acyclic_algorithm(B2,ω,h1,γ1)
@info "Sparse matrix and uniform coupling: OK!"

xxx3 = run_acyclic_algorithm(B1,ω,h2,γ2)
@info "Dense matrix and heterogeneous couplings: OK!"

xxx4 = run_acyclic_algorithm(B2,ω,h2,γ2)
@info "Sparse matrix and heterogeneous couplings: OK!"





