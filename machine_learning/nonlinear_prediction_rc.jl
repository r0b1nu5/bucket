using DelimitedFiles, PyPlot

include("lorenz.jl")
include("henon-heiles.jl")
include("reservoir.jl")

model = "lorenz"
new_reservoir = false

if model == "lorenz"
# #=
#	xs1 = readdlm("data1/xs1.csv",',')
#	xs2 = readdlm("data1/xs2.csv",',')
#	xs3 = readdlm("data1/xs3.csv",',')
	xs1 = readdlm("data1/xs4.csv",',')
# =#
	n = 3
elseif model == "hh"
#	xs1 = readdlm("data1/hh_xys0_chaos.csv",',')
#	xs1 = readdlm("data1/hh_xys0_periodic.csv",',')
#	xs1 = readdlm("data1/hh_xys1_periodic.csv",',')
	
	xs1 = xys2

	n = 4
end

#=
T0 = 1 # starting time
Tt = 8001 # training time
=#
Tp = 1000 # time between training and prediction
DT = 1000 # prediction time

thr = 1e-1

N = 1000
m = 10*N
rho = 1.5
dT = 1

sig = 1.

a = 1.
xi = 1.
beta = .01

if new_reservoir
	A = A_gen(N,m,rho)
	Win = Win_gen(n,N,sig)
end

Wout,c,r = reservoir_training((xs1[:,T0:T0+Tt-1],xs1[:,T0+1:T0+Tt]),A,Win,a,xi,dT,beta)

if DT > 0
	rr = reservoir_tanh(r[:,end],xs1[:,T0+Tt:T0+Tt+DT],A,Win,a,xi)
else
	rr = r
end

ss = Array{Float64,2}(undef,length(xs1[:,T0+Tt]),0)
ss = [ss xs1[:,T0+Tt+DT]]
rs = Array{Float64,2}(undef,length(r[:,end]),0)
rs = [rs rr[:,end]]

for t in 1:Tp
	global ss,rs
	s,r = reservoir_prediction(ss[:,[size(ss)[2],]],Wout,c,rs[:,end],A,Win,a,xi)
	ss = [ss s]
	rs = [rs r]
end

mas = maximum(xs1[:,(T0+Tt+DT):(T0+Tt+DT+Tp)],dims=2) 
mis = minimum(xs1[:,(T0+Tt+DT):(T0+Tt+DT+Tp)],dims=2)
amps = mas - mis
diff = xs1[:,(T0+Tt+DT):(T0+Tt+DT+Tp)] - ss
errs = abs.(diff)./repeat(amps,1,Tp+1)
merr = [maximum(errs[:,1:i]) for i in 1:Tp+1]

Tb = minimum(setdiff((1:Tp+1).*(merr .> thr),[0,]))

 #=
fignum = rand(1000:9999)
for i in 1:n
	figure(fignum)
	subplot(n+1,1,i)
	PyPlot.plot((T0+Tt+DT):(T0+Tt+DT+Tp),xs1[i,(T0+Tt+DT):(T0+Tt+DT+Tp)])
	PyPlot.plot((T0+Tt+DT):(T0+Tt+DT+Tp),ss[i,:],"--")
	PyPlot.plot([T0+Tt+DT+Tb,T0+Tt+DT+Tb],[minimum(xs1[i,(T0+Tt+DT):(T0+Tt+DT+Tp)]),maximum(xs1[i,(T0+Tt+DT):(T0+Tt+DT+Tp)])],"--k")
	axis([T0+Tt+DT-100,T0+Tt+DT+Tp+100,mis[i],mas[i]])
	subplot(n+1,1,n+1)
	PyPlot.semilogy((T0+Tt+DT):(T0+Tt+DT+Tp),errs[i,:],color="C$(i+1)")
end
figure(fignum)
subplot(n+1,1,n+1)
PyPlot.semilogy((T0+Tt+DT):(T0+Tt+DT+Tp),merr,"--k")
PyPlot.plot([T0+Tt+DT+Tb,T0+Tt+DT+Tb],[1e-3,1.1],"--k")
axis([T0+Tt+DT-100,T0+Tt+DT+Tp+100,1e-3,1.1])
# =#




