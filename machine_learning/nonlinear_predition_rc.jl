using DelimitedFiles, PyPlot

include("lorenz.jl")
include("reservoir.jl")

model = "lorenz"
new_reservoir = false

if model == "lorenz"
# #=
	xs1 = readdlm("data1/xs1.csv",',')
	xs2 = readdlm("data1/xs2.csv",',')
	xs3 = readdlm("data1/xs3.csv",',')
# =#
	n = 3
end

T = 4001
Tp = 200
DT = 0

if new_reservoir
	N = 1000
	m = 10*N
	rho = 1.
	A = A_gen(N,m,rho)
	dT = 1

	sig = 1.
	Win = Win_gen(n,N,sig)

	a = 1.
	xi = 1.
	beta = .01
end

Wout,c,r = reservoir_training((xs1[:,1:T-1],xs1[:,2:T]),A,Win,a,xi,dT,beta)

if T > 0
	rr = reservoir_tanh(r[:,end],xs1[:,T:T+DT],A,Win,a,xi)
else
	rr = r
end

ss = Array{Float64,2}(undef,length(xs1[:,T]),0)
ss = [ss xs1[:,T+DT]]
rs = Array{Float64,2}(undef,length(r[:,end]),0)
rs = [rs rr[:,end]]

for t in 1:Tp
	global ss,rs
	s,r = reservoir_prediction(ss[:,[size(ss)[2],]],Wout,c,rs[:,end],A,Win,a,xi)
	ss = [ss s]
	rs = [rs r]
end

figure()
subplot(3,1,1)
PyPlot.plot(T+DT:(T+DT+Tp),xs1[1,T+DT:(T+DT+Tp)])
PyPlot.plot(T+DT:(T+DT+Tp),ss[1,:],"--")
subplot(3,1,2)
PyPlot.plot(T+DT:(T+DT+Tp),xs1[2,T+DT:(T+DT+Tp)])
PyPlot.plot(T+DT:(T+DT+Tp),ss[2,:],"--")
subplot(3,1,3)
PyPlot.plot(T+DT:(T+DT+Tp),xs1[3,T+DT:(T+DT+Tp)])
PyPlot.plot(T+DT:(T+DT+Tp),ss[3,:],"--")


