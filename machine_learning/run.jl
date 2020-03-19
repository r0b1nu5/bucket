include("reservoir.jl")

xs = readdlm("data1/xs4.csv",',')

Tp = 1000 # predition time
DT = 1000 # time between training and prediction

N = 2000 # reservoir size

n = 3
m = 10*N
rho = 1.5
sig = 1.
a = 1.
xi = 1.
dT = 1
beta = .01

thr = .1

A = A_gen(N,m,rho)
Win = Win_gen(n,N,sig)

# Compute Wout and breaktime with respect to training time.
# #=
Tts = Array(101:100:4001) # training times

for i in 1:length(Tts)
	global Tt = Tts[i]
	global T0 = 8002 - Tt
	Wout,c,r = reservoir_training((xs[:,T0:T0+Tt-1],xs[:,T0+1:T0+Tt]),A,Win,a,xi,dT,beta)
	writedlm("data1/Wout$(i)_vs_Tt_N$N.csv",Wout,',')

	if DT > 0
		rr = reservoir_tanh(r[:,end],xs[:,T0+Tt:T0+Tt+DT],A,Win,a,xi)
	else
		rr = r
	end

	ss,rs = reservoir_prediction_self(xs[:,T0+Tt+DT],rr[:,end],Tp,Wout,c,A,Win,a,xi)


	Tb = breaktime(thr,xs[:,T0+Tt+DT:end],ss)
	writedlm("data1/Tb$(i)_vs_Tt_N$N.csv",Tb,',')

	@info "Tt = $Tt, i = $i done."
end
# =#

Tt = 1001
T0 = 8002 - Tt
Ns = Array(200:200:4000)

# #=
for i in 1:length(Ns)
for k in 1:50
	global N = Ns[i]
	m = 10*N
	
	global A = A_gen(N,m,rho)
	global Win = Win_gen(n,N,sig)

	Wout,c,r = reservoir_training((xs[:,T0:T0+Tt-1],xs[:,T0+1:T0+Tt]),A,Win,a,xi,dT,beta)
	writedlm("data1/Wout$(k).$(i)_vs_N_Tt$(Tt).csv",Wout,',')

	if DT > 0
		rr = reservoir_tanh(r[:,end],xs[:,T0+Tt:T0+Tt+DT],A,Win,a,xi)
	else
		rr = r
	end

	ss,rs = reservoir_prediction_self(xs[:,T0+Tt+DT],rr[:,end],Tp,Wout,c,A,Win,a,xi)

	Tb = breaktime(thr,xs[:,T0+Tt+DT:end],ss)
	writedlm("data1/Tb$(k).$(i)_vs_N_Tt$(Tt).csv",Tb,',')
	
	@info "N = $N, i = $i done."
end
end
# =#



