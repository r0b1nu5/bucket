using PyPlot, Dates

include("tools.jl")

t0 = time()

v = [0,0,0,0,0,0,0,1.]
col = "C7"

ρs = LinRange(-2,2,200)
N = 100000

L = [3 -1 -1 -1;-1 2 -1 0;-1 -1 2 0;-1 0 0 1.]
L = readdlm("ntw_data/ntw10_L.csv",',')

b,w = L2B(L)
B = [b -b]
Bout = B.*(B .> 0.)

n,m = size(B)

Id = diagm(0 => ones(m))

λs = Array{Float64,2}(undef,m,0)

PD = Id - B'*pinv(Bout*B')*Bout
Q = 2*PD'*PD + Id - PD - PD'
#λs = [λs eigvals(Q)]
#ρ0 = λs[m,1]/λs[1,1]
#ρ = λs[m,1]/λs[1,1]
vmin = ones(m)
imin = 1

αs = LinRange(0,10,500)
for α in αs
	PD = Id - B'*pinv(Bout*B')*Bout

	Q = (1+α)*PD'*PD + α*(Id - PD - PD')
	
	λ = real.(eigvals(Q))
	global λs = [λs λ]


end

PyPlot.plot(αs[2:end],λs[m,2:end]./abs.(λs[1,2:end]))

@info "$(time() - t0)''"

