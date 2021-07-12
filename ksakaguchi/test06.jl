using PyPlot

include("iterations.jl")

L = readdlm("ntw_data/ntw10_L.csv",',')
b,w = L2B(L)
bd = pinv(b)
B = [b -b]
Bout = B.*(B .> 0)
n,m = size(B)
m2 = Int(m/2)
P = cycle_proj(b,ones(m2))
Id = diagm(0 => ones(m2))

mi = .0
ma = 1.
de = ma - mi

n_iter = 1000

for i in 1:n_iter
	d1 = de*rand(m2) .+ mi
	d2 = d1
	d2 += .06*rand(m2) .- .03
	d2 .-= min(minimum(d2),0.)
#	d2 = de*rand(m2) .+ mi
	d0 = (d1 + d2)/2
	dm = (d1 - d2)/2
#	d0 = mean([d1;d2])
#	dm = [d1;d2] .- d0
	
	D = [diagm(0 => d1); diagm(0 => -d2)]
	ρ1 = .1
	ρ2 = .1
	
	M = Id - ρ1*bd*Bout*D - ρ2*P
	
	λs = eigvals((M+M')/2)
	
	if maximum(λs) > 1.
		c = "C3"
		global R = [R d1 d2]
	else
		c = "C0"
	end
		
	PyPlot.plot(maximum(λs),maximum(abs.(dm)./d0),".",color=c)
end

