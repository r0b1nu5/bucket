using DelimitedFiles

include("acyclic_algorithm_repo.jl")

h = ((x,a) -> (sin(x - a) + sin(a)))

ω = .6*rand(8)
ω .- mean(ω)

L1 = readdlm("ntw_data/ntw8_acyc_L.csv",',')

IJV = readdlm("ntw_data/ntw8_acyc_Lsp.csv",',')
L2 = sparse(Int64.(IJV[:,1]),Int64.(IJV[:,2]),IJV[:,3])

n = 8
m = 7

α = .3
γ1 = π/2 - α
H1 = (f -> (-sin(asin(f - sin(α)) + 2*α) + sin(α)))
hγ1 = (max(h(-γ1,α),-h(γ1,α)),min(h(γ1,α),-h(-γ1,α)))

αs = .4*rand(14) .+ .1
γ2 = π/2 .- αs
H2 = Vector{Function}()
hγ2 = Vector{Tuple{Float64,Float64}}()
for i in 1:m
	push!(H2,(f -> (-sin(asin(f - sin(αs[i+m])) + αs[i] + αs[i+m]) + sin(αs[i]))))
	push!(hγ2,(max(h(-γ2[i],αs[i]),-h(γ2[i+m],αs[i+m])),min(h(γ2[i],αs[i]),-h(-γ2[i+m],αs[i+m]))))
end
for i in m+1:2*m
	push!(H2,(f -> (-sin(asin(f - sin(αs[i-m])) + αs[i] + αs[i-m]) + sin(αs[i]))))
	push!(hγ2,(max(h(-γ2[i],αs[i]),-h(γ2[i-m],αs[i-m])),min(h(γ2[i],αs[i]),-h(-γ2[i-m],αs[i-m]))))
end


xxx1 = run_acyclic_algorithm(L1,H1,ω,hγ1)
@info "Dense matrix and uniform coupling: OK!"

xxx2 = run_acyclic_algorithm(L2,H1,ω,hγ1)
@info "Sparse matrix and uniform coupling: OK!"

xxx3 = run_acyclic_algorithm(L1,H2,ω,hγ2)
@info "Dense matrix and heterogeneous couplings: OK!"

xxx4 = run_acyclic_algorithm(L2,H2,ω,hγ2)
@info "Sparse matrix and heterogeneous couplings: OK!"





