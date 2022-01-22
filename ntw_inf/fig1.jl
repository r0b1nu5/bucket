using PyPlot, DelimitedFiles, Statistics, SparseArrays

include("line_dist.jl")

ntw = "pegase1354"
ρ = .1

Lsp = readdlm("ntw_data/"*ntw*"_lap_mat_sp.csv",',')
L = sparse(Int64.(Lsp[:,1]),Int64.(Lsp[:,2]),Lsp[:,3])
B,w,Bt = L2B(L)
n,m = size(B)
Ld = pinv(Matrix(L))

P = ρ*rand(n)
P .-= mean(P)

l = 6

θ0s,it = kuramoto_sine(ntw,L,P,Ld*P,1,0.,.1,0.,false,20000,1e-5,.1)
θ0 = θ0s[:,end]

J = spzeros(n,n)
for k in 1:size(Lsp)[1]
	i = Int64(Lsp[k,1])
	j = Int64(Lsp[k,2])
	if i < j
		J[i,j] = -Lsp[k,3]*cos(θ0[i]-θ0[j])
		J[j,i] = -Lsp[k,3]*cos(θ0[i]-θ0[j])
	end
end
J = J - spdiagm(0 => vec(sum(J,dims=2)))
λs = eigvals(Matrix(J))
λ2 = -λs[end-1]
λn = -λs[1]

τ = 1/(10*λn)

θs,it = kuramoto_sine(ntw,L,P,θ0,l,1.,1/τ,0.,true,1000,-1e-5,min(.01,τ/10))
ψs = L*θs

N88 = [92,93,109,239,216,527]

figure()
for i in [88;N88;150:200]
	PyPlot.plot((501:1000)*τ,θs[i,501:1000] .- mean(θs[i,501:1000]))
end
xlabel("t")
ylabel("x")

figure()
for i in [88;N88;150:200]
	PyPlot.plot((501:1000)*τ,ψs[i,501:1000] .- mean(ψs[i,501:1000]))
end
xlabel("t")
ylabel("ψ")

