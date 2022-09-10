using PyPlot, DelimitedFiles

obj = vec(readdlm("data/ntw20_obj_vs_freq.csv",','))

ks = 1:150
N = 50000
τ = .01
T = N*τ

PyPlot.plot(ks/T,obj)
xlabel("freq")
ylabel("obj")



