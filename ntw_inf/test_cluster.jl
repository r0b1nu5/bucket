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

function test_cluster(L,n,P,th0,dP0,tau0,h,t,T,dm)
	ti = time()

	x = noise_inf(L,P,th0,dP0,tau0,t,T,dm,1e-8,h,true)

	return time() - ti
end

function test_cluster(L,n,P,th0,dP0,tau0,h,t,T,dm)
	ti = time()

	x = noise_inf(L,P,th0,dP0,tau0,t,T,dm,1e-8,h,true)

	return time() - ti
end


