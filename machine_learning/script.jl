include("reservoir.jl")

a = 1.
xi = 1.
dT = 1
beta = .01
rho = 1.5
N = 2000
m = 10*N
n = 3
sig = 1.

Tbsss = Array{Array{Int64,2},1}()
Xs = Array{Array{Float64,2},1}()
for i in 2:4
	push!(Xs,readdlm("data1/xs$i.csv",','))
end

for i in 1:100
	Tbss = Array{Int64,2}(undef,3,0)
	A = A_gen(N,m,rho)
	Win = Win_gen(n,N,sig)
	for xs in Xs
		Wout,c,r = reservoir_training((xs[:,T0:T0+Tt-1],xs[:,T0+1:T0+Tt]),A,Win,a,xi,dT,beta)
		rr = reservoir_tanh(r[:,end],xs[:,T0+Tt:T0+Tt+DT],A,Win,a,xi)
		ss,rs = reservoir_prediction_self(xs[:,T0+Tt+DT],rr[:,end],Tp,Wout,c,A,Win,a,xi)
		Tbss = [Tbss [breaktime(.1,xs[:,T0+Tt+DT:end],ss),breaktime(.05,xs[:,T0+Tt+DT:end],ss),breaktime(.02,xs[:,T0+Tt+DT:end],ss)]]
	end
	push!(Tbsss,Tbss)
end

mTbs = sum(Tbsss)/100


