using DelimitedFiles

n = 120 

Asp = readdlm("data/uk_adj_mat.csv",',')
Lsp = readdlm("data/lap_mat.csv",',')
L = zeros(n,n)

for i in 1:size(Lsp)[1]
	L[Int(Lsp[i,1]),Int(Lsp[i,2])] = Lsp[i,3]
end

m = ones(n)
d = ones(n)



