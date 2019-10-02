using PyPlot, DelimitedFiles, Statistics, LinearAlgebra

ntw = "ntw20"
n = 20
sources = Array(1:n)
a0 = .2
f0 = .00945
p0 = pi/10

for s in sources
	afp = readdlm("data/"*ntw*"_afp_$s.csv",',')
	a = afp[:,1]
	f = afp[:,2]
	p = afp[:,3]

	AA = sortslices([a 1:n],dims=1,rev=true)
	
	figure(50)
	for i in 1:n
		if Int(AA[i,2]) == s
			PyPlot.plot(s,log(AA[i,1]/AA[1,1]),"x",color="C2")
		else
			PyPlot.plot(s,log(AA[i,1]/AA[1,1]),"x",color="C3")
		end
	end

	figure(51)
	PyPlot.semilogy(s,abs(a[Int(AA[1,2])] - a0)/a0,"x",color="C0")
	PyPlot.semilogy(s,abs(f[Int(AA[1,2])] - f0)/f0,"x",color="C1")
	PyPlot.semilogy(s,abs(p[Int(AA[1,2])] - p0)/p0,"x",color="C4")

end

figure(51)
legend(["a","f","φ"])








