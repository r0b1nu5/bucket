using PyPlot, DelimitedFiles

σs = [2.825,2.9,2.975,3.05,3.125,3.25]

figure(figsize=(15,6))

c = 0
for σ in σs
	global c += 1
	subplot(2,3,c)
	x = readdlm("data-pj/sigma_$σ.dat")
	PyPlot.plot(x[1:20000,1],x[1:20000,2])
	title("σ = $σ")
	xlabel("t [a.u.]")
	ylabel("N")
end



