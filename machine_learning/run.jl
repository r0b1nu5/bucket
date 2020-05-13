include("reservoir.jl")

xs_file = "_pj3"
writedlm("data1/last_run_xs_file.csv",xs_file,',')

xs = readdlm("data1/xs"*xs_file*".csv",',')

Tp = 1000 # predcition time
DT = 1000 # time between training and prediction

Tti = 2001
Ttf = 2001
dTt = 50
Ni = 2000
Nf = 2000
dN = 200
k0 = 20
ks = 80

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

Tts = Array(Tti:dTt:Ttf)
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

 #=
# Compute breaktime with respect to reservoir size.

Ns = Array(Ni:dN:Nf)
Tt0 = 2001
writedlm("data1/last_run_Ns.csv",Ns,',')
writedlm("data1/last_run_Tt0.csv",Tt0,',')

Tt = Tt0
T0 = 8002 - Tt

for i in 1:length(Ns)
for k in 1:ks
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

# #=
# Compute breaktime with respect to reservoir size and training time

Ns = Array(Ni:dN:Nf)
Tts = Array(Tti:dTt:Ttf)
writedlm("data1/last_run_Ns.csv",Ns,',')
writedlm("data1/last_run_Tts.csv",Tts,',')

for i in 1:length(Ns)
	for k in (k0+1):(k0+ks)
		global N = Ns[i]
		m = 10*N

		global A = A_gen(N,m,rho)
		global Win = Win_gen(n,N,sig)

		for j in 1:length(Tts)
			global Tt = Tts[i]
			global T0 = 8002 - Tt

			Wout,c,r = reservoir_training((xs[:,T0:T0+Tt-1],xs[:,T0+1:T0+Tt]),A,Win,a,xi,dT,beta)
			writedlm("data1/Wout"*xs_file*"_$(k).$(i).$(j)_vs_N_vs_Tt.csv",Wout,',')

			if DT > 0
				rr = reservoir_tanh(r[:,end],xs[:,T0+Tt:T0+Tt+DT],A,Win,a,xi)
			else
				rr = r
			end

			ss,rs = reservoir_prediction_self(xs[:,T0+Tt+DT],rr[:,end],Tp,Wout,c,A,Win,a,xi)

			for thr in thrs
				Tb = breaktime(thr,xs[:,T0+Tt+DT:end],ss)
				writedlm("data1/Tb"*xs_file*"_$(k).$(i).$(j)_vs_N_vs_Tt_thr$(thr).csv",Tb,',')
			end

			@info "N = $N, Tt = $Tt, k = $k done."
		end
	end
end

# =#


