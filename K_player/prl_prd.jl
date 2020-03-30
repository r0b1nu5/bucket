using DelimitedFiles, PyPlot

dat_prl = readdlm("data/prl.txt",'\t')
dat_prd = readdlm("data/prd.txt",'\t')

mi_prl = 63
ma_prl = 102
mi_prd = 72
ma_prd = 111

prl = Dict{Int64,Array{String,1}}()
prd = Dict{Int64,Array{String,1}}()
prl_prd = Dict{Tuple{Int64,Int64},Array{String,1}}()
n_prl_prd = Dict{Tuple{Int64,Int64},Int64}()
matri = zeros(ma_prl-mi_prl+1,ma_prd-mi_prd+1)

for i in mi_prl:ma_prl
	temp = Array{String,1}()
	for j in 1:size(dat_prl)[1]
		if dat_prl[j,2] == i
			push!(temp,dat_prl[j,1])
		end
	end
	prl[i] = temp
end

for i in mi_prd:ma_prd
	temp = Array{String,1}()
	for j in 1:size(dat_prd)[1]
		if dat_prd[j,2] == i
			push!(temp,dat_prd[j,1])
		end
	end
	prd[i] = temp
end

for i in keys(prl)
	for j in keys(prd)
		x = intersect(prl[i],prd[j])
		l = length(x)
		if l > 0
			prl_prd[(i,j)] = x
			n_prl_prd[(i,j)] = l
			matri[i-mi_prl+1,j-mi_prd+1] = l
		end
	end
end

matri = permutedims(matri)
ccmap = "plasma"

matshow(matri,cmap=ccmap)
colorbar()
xticks((0:5:39),(mi_prl:5:ma_prl))
xlabel("PRL: # of articles")
yticks((0:5:37),(mi_prd:5:ma_prd))
ylabel("PRD: # of articles")
title("# of authors")

matshow(max.(zeros(size(matri)),log10.(matri.+.25)),cmap=ccmap)
#matshow(log10.(matri),cmap=ccmap)
colorbar()
xticks((2:5:39),(mi_prl+2:5:ma_prl))
xlabel("PRL: # of articles")
yticks((3:5:38),(mi_prd+3:5:ma_prd))
ylabel("PRD: # of articles")
title("log(# of authors)")

println("(96,77): CMS experiment at LHC, CERN")
println("(66,104): ATLAS experiment at LHC, CERN")


