using PyPlot, DelimitedFiles

ids = [2,3,4,5,6,7]

res = Dict{Int64,Matrix{Float64}}(i => readdlm("eeg-data/relative-error-x$(i)x.csv",',') for i in ids)

ss = ["001","002"]
rs = ["01","02"]

figure("Relative error",(6,3))

for i in 1:length(ss)
	for j in 1:length(rs)
		re = [res[id][j,i] for id in ids]
		PyPlot.plot(ids,re,"o-",label="S"*ss[i]*"R"*rs[j])
	end
end

xlabel("order")
ylabel("relative error")
legend()

