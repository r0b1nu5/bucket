using DelimitedFiles

idx = 45
freq = 100.

n = 120 

Lsp = readdlm("data/lap_mat.csv",',')
L = zeros(n,n)

for i in 1:size(Lsp)[1]
	L[Int(Lsp[i,1]),Int(Lsp[i,2])] = Lsp[i,3]
end

m = ones(n)
d = ones(n)

c = zeros(n)
c[idx] = 1.
f = zeros(n)
f[idx] = freq
phi = zeros(n)
phi[idx] = pi/10
forc = (c,f,phi)

sig = ones(n)

ntw = "uk$idx"


