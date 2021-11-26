using DelimitedFiles

include("acyclic_algorithm_repo.jl")

ω = .6*rand(8)
ω .- mean(ω)
ω = zeros(8)

L1 = readdlm("ntw_data/ntw8_acyc_L.csv",',')

IJV = readdlm("ntw_data/ntw8_acyc_Lsp.csv",',')
L2 = sparse(Int64.(IJV[:,1]),Int64.(IJV[:,2]),IJV[:,3])

α = .3
γ1 = π/2 - α
H1 = (x -> (sin(x - α) - sin(α)))
hγ1 = (H1(-γ1),H1(γ1))

αs = .4*rand(14) .+ .1
γ2 = π/2 .- αs
H2 = Vector{Function}()
hγ2 = Vector{Tuple{Float64,Float64}}()
for i in 1:length(αs)
	push!(H2,(x -> (sin(x - αs[i]) - sin(αs[i]))))
	push!(hγ2,(H2[i](-γ2[i]),H2[i](γ2[i])))
end


xxx1 = run_acyclic_algorithm(L1,H1,ω,hγ1)
@info "Dense matrix and uniform coupling: OK!"

xxx2 = run_acyclic_algorithm(L2,H1,ω,hγ1)
@info "Sparse matrix and uniform coupling: OK!"

xxx3 = run_acyclic_algorithm(L1,H2,ω,hγ2)
@info "Dense matrix and heterogeneous couplings: OK!"

xxx4 = run_acyclic_algorithm(L2,H2,ω,hγ2)
@info "Sparse matrix and heterogeneous couplings: OK!"





