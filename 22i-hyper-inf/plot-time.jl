using DelimitedFiles, PyPlot

ns = [16,32,64,128]
runs = ["995","996","997","998","999"]
Ts = [1000,2000]

figure(1,(4,3.5))
for T in Ts
	tt = Float64[]
	for n in ns
		t = Float64[]
		for run in runs
			push!(t,readdlm("data/kuramoto-"*run*"-n$n-time$T-1.csv",',')[1])
		end
		push!(tt,mean(t))
	end
	figure(1)
	PyPlot.loglog(ns,tt,"o",label="T = $T")
end
xlabel("N")
ylabel("computation time [s]")
legend()
axis([9,200,.01,1000])

