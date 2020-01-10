using Distributed, DelimitedFiles

@everywhere include("load_ntw.jl")
@everywhere include("kuramoto.jl")

n_runs = 5
n_sruns = 1000
ntw = "ntw3"
ee = -.00005

params = Array{Tuple{Int64,Int64,String,Float64},1}()
for i in 1:n_runs
	push!(params,(i,n_sruns,ntw,ee))
end

@everywhere function run(param::Tuple{Int64,Int64,String,Float64})
	r_id,n_sruns,ntw,ee = param

	L,P = load_ntw(ntw,ee)

	n = length(P)
	
	its = Array{Int64,1}()
	th0s = Array{Float64,2}(undef,n,0)

	for i in 1:n_sruns
		@info "$r_id, $i"

		th0 = [2*pi*rand(n-1);0]
		t,d,it = kuramoto(L,P,th0,false,false)
		push!(its,it)
		th0s = [th0s th0]
	end

	writedlm("data/"*ntw*"_$(r_id)_th0s.csv",th0s,',')
	writedlm("data/"*ntw*"_$(r_id)_its.csv",its,',')
end

pmap(run,params)



	



