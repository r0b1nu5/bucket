include("noisy_inf.jl")

Lsp = readdlm("ntws_data/ntw20_lap_mat_sp.csv",',')
L = sparse(Lsp[:,1],Lsp[:,2],Lsp[:,3])

n = size(L)[1]
P = zeros(n)
th0 = zeros(n)
dP0 = .1
tau0 = .001
h = .0001
t = 1
T = 300000
dm = 1

dPs = readdlm("data1/dPs.csv",',')

function test_cluster(L,n,P,th0,dP0,dPs,h)
	ti = time()

	x = noise_inf_prerand(L,P,th0,dP0,dPs,1e-8,h,true)

	return time() - ti
end



