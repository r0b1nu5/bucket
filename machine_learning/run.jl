include("reservoir.jl")

xs_file = 1
writedlm("data1/last_run_xs_file.csv",xs_file,',')

xs = readdlm("data1/xs$(xs_file).csv",',')

Tp = 1000 # predition time
DT = 1000 # time between training and prediction

n = 3
rho = 1.5
sig = 1.
a = 1.
xi = 1.
dT = 1
beta = .01

thrs = [.1,.05,.02]


 #=
# Compute Wout and breaktime with respect to training time.

Tts = Array(51:50:2001)
N = 2000
writedlm("data1/last_run_Tts.csv",Tts,',')
writedlm("data1/last_run_N0.csv",N,',')

m = 10*N
A = A_gen(N,m,rho)
Win = Win_gen(n,N,sig)

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

	for thr in thrs
		Tb = breaktime(thr,xs[:,T0+Tt+DT:end],ss)
		writedlm("data1/Tb$(i)_vs_Tt_N$(N)_thr$(thr).csv",Tb,',')
	end

	@info "Tt = $Tt, i = $i done."
end
# =#

# #=
# Compute breaktime with respect to reservoir size.

Ns = Array(200:200:2000)
Tt0 = 1001
writedlm("data1/last_run_Ns.csv",Ns,',')
writedlm("data1/last_run_Tt0.csv",Tt0,',')

Tt = Tt0
T0 = 8002 - Tt

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

	for thr in thrs
		Tb = breaktime(thr,xs[:,T0+Tt+DT:end],ss)
		writedlm("data1/Tb$(k).$(i)_vs_N_Tt$(Tt)_thr$(thr).csv",Tb,',')
	end
	
	@info "N = $N, i = $i done."
end
end
# =#


