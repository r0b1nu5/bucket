include("kuramoto.jl")
include("transfo.jl")

L = readdlm("ntws_data/ntw12_lap_mat.csv",',')

n = 12

P = 1.5*rand(n)
P .-= mean(P)

a0 = .2
w0 = .01
p0 = 0.

T = 20000
h = .1

Ic = [10,11,12]

for ij in [(7,8),(2,10),(11,12)]
	ths,dths = kuramoto_transfo(L,P,zeros(n),ij,(a0,w0,p0),T,true,1e-6,h)
	ps = psi(L,Ic,ths)

	figure("(i,j) = ($(ij[1]),$(ij[2]))")
	for i in 1:9
		PyPlot.plot(h*(1:T),ps[i,:] .- ps[i,1])
	end
end

for i in [6,11]
	a = zeros(n)
	a[i] = a0
	w = zeros(n)
	w[i] = w0
	p = zeros(n)

	ths = kuramoto_sine(L,P,zeros(n),zeros(n),zeros(n),zeros(n),1)

	th1 = ths[:,end]

	ths = kuramoto_sine(L,P,th1,a,w,p,T,T)
	ps = psi(L,Ic,ths)

	figure("i = $i")
	for i in 1:9
		PyPlot.plot(h*(1:T),ps[i,1:end-1] .- ps[i,1])
	end
end




