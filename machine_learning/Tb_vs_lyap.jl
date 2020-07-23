using PyPlot, DelimitedFiles

include("reservoir.jl")
include("lyap.jl")

files = ["hh_xys01", "hh_xys02", "hh_xys03", "hh_xys04"]

for fi in files
	xs = readdlm("data1/"*fi*".csv",',')

	T0 = 2000
	Tt = 8000
	DT = 1000
	Tp = 2000
	N = 2000
	k0 = 100
	ks = 20

	n,T = size(xs)
	rho = 1.5
	sig = 1.
	a = 1.
	xi = 1.
	dT = 1
	beta = .01
	m = 10*N

	thrs = [.1,.05,.02]

	l = readdlm("data1/"*fi*"_lyap.csv",',')[1]

	for k in k0+1:k0+ks
		A = A_gen(N,m,rho)
		Win = Win_gen(n,N,sig)

		Wout,c,r = reservoir_training((xs[:,T0:T0+Tt-1],xs[:,T0+1:T0+Tt]),A,Win,a,xi,dT,beta)

		rr = reservoir_tanh(r[:,end],xs[:,T0+Tt:T0+Tt+DT],A,Win,a,xi)

		ss,rs = reservoir_prediction_self(xs[:,T0+Tt+DT],rr[:,end],Tp,Wout,c,A,Win,a,xi)

		for thr in thrs
			Tb = breaktime(thr,xs[:,T0+Tt+DT:end],ss)
			writedlm("data1/lyap_Tb_"*fi*"_k$(k)_thr$thr.csv",[l,Tb],',')
		end
	end
end






