using DelimitedFiles, PyPlot, Statistics

include("../L2B.jl")
include("kuramoto.jl")
include("tools.jl")

L = readdlm("ntw_data/ntw10_L.csv",',')
B,w = L2B(L)
n,m = size(B)

ρ = .2

P = ρ*rand(n)
P .-= mean(P)

niter = 100

θis = zeros(n,0)

for i in 1:niter
	global n,L,P,θis

	θ0 = 2π*rand(n)

	θs,xxx = kuramoto(L,P,θ0,true)
	n,T = size(θs)

	λs = zeros(n,0)
	for t in 1:T
		λs = [λs eigvals(jacobian(θs[:,t],L))]
	end

	sig = vec(sum(λs .> 1e-10,dims=1))

	dsig = sig[2:end] - sig[1:end-1]

	if maximum(dsig) > .1
		θis = [θis θ0]
		@info "BINGO! #$(size(θis)[2])"
	end
end





