using DelimitedFiles, PyPlot

include("lorenz.jl")
include("henon-heiles.jl")
include("reservoir.jl")

model = "hh"
new_reservoir = true

if model == "lorenz"
# #=
	xs1 = readdlm("data1/xs1.csv",',')
	xs2 = readdlm("data1/xs2.csv",',')
	xs3 = readdlm("data1/xs3.csv",',')
# =#
	n = 3
elseif model == "hh"
#	xs1 = readdlm("data1/hh_xys0_chaos.csv",',')
#	xs1 = readdlm("data1/hh_xys0_periodic.csv",',')
#	xs1 = readdlm("data1/hh_xys1_periodic.csv",',')
	
	xs1 = xys2

	n = 4
end

T = 4001
Tp = 1000
DT = 1000

if new_reservoir
	N = 1000
	m = 10*N
	rho = 1.5
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
for i in 1:n
	subplot(n,1,i)
	PyPlot.plot(T+DT:(T+DT+Tp),xs1[i,T+DT:(T+DT+Tp)])
	PyPlot.plot(T+DT:(T+DT+Tp),ss[i,:],"--")
end


