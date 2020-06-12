using Distributed

@everywhere include("reservoir.jl")

xs_file = "_pj3"
xs = readdlm("data1/xs"*xs_file*".csv",',')

Tp = 1000 # prediction time
DT = 1000 # time between training and prediction

id = rand(10000:99999) # Run's ID

Tti = 101
Ttf = 2001
dTt = 200
Ni = 200
Nf = 2000
dN = 200
k0 = 0
ks = 100

n = 3
rho = 1.5
sig = 1.
a = 1.
xi = 1.
dT = 1
beta = .01

thrs = [.1,.05,.02]

writedlm("data1/run_$(id)_params.jl",[Tp,DT,Tti,Ttf,dTt,Ni,Nf,dN,k0,ks,n,rho,sig,a,xi,dT,beta,thrs],',')

Ns  = Array(Ni:dN:Nf)
Tts = Array(Tti:dTt:Ttf)

args = Array{Tuple{Int64,Int64,Int64,Int64,Array{Float64,2},Array{Float64,2},Array{Float64,2},Float64,Float64,Int64,Float64,Int64,Int64,Array{Float64,1}},1}()

for i in 1:length(Ns)
	global args
	for k in (k0+1):(k0+ks)
		N = Ns[i]
		m = 10*N
		A = A_gen(N,m,rho)
		Win = Win_gen(n,N,sig)
		for j in 1:length(Tts)
			Tt = Tts[j]
			T0 = 8002 - Tt
			push!(args,(id,i,j,k,xs[:,T0:T0+Tt+DT+Tp],A,Win,a,xi,dT,beta,DT,Tp,thrs))
		end
	end
end

@everywhere function run_script(arg::Tuple{Int64,Int64,Int64,Int64,Array{Float64,2},Array{Float64,2},Array{Float64,2},Float64,Float64,Int64,Float64,Int64,Int64,Array{Float64,1}})
	id,i,j,k,xs,A,Win,a,xi,dT,beta,DT,Tp,thrs = arg
	
	n,T = size(xs)
	N = size(Win)[1]
	Tt = T - DT - Tp

	Wout,c,r = reservoir_training((xs[:,1:T-Tp-DT-1],xs[:,2:T-Tp-DT]),A,Win,a,xi,dT,beta)

	subid = (i,j,k)
	if DT > 0
		rr = reservoir_tanh(r[:,end],xs[:,T-Tp-DT:T-Tp],A,Win,a,xi,subid)
	else
		rr = r
	end

	ss,rs = reservoir_prediction_self(xs[:,T-Tp],rr[:,end],Tp,Wout,c,A,Win,a,xi,subid)

	for thr in thrs
		Tb = breaktime(thr,xs[:,T-Tp:end],ss)
		writedlm("data1/Tb_$(k).$(i).$(j)_$(thr)_$(id).csv",Tb,',')
	end

	@info "Run $(id): N = $N, Tt = $Tt, k = $k done."
end

pmap(run_script,args)


