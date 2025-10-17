using DelimitedFiles, PowerModels, PyPlot

include("tools.jl")

function gen_2p_correlated_time_series(i1::Int64, i2::Int64, r::Float64=0., l::Int64=86400)
	# Loading network data
	network_data = load_testcase_1()

	# Generating power time series
	μp = vec(readdlm("rv/data-marc/mup.csv",','))
	σp = vec(readdlm("rv/data-marc/sigp.csv",','))
	μq = vec(readdlm("rv/data-marc/muq.csv",','))
	σq = vec(readdlm("rv/data-marc/sigq.csv",','))
	times = vec(readdlm("rv/data-marc/times.csv",','))[1:l]

	head = ["Times" (1:55)']

	n = length(μp)
	R = zeros(1,n); R[1,i1] = R[1,i2] = r

	x0 = randn(l) # Reference time series
	p = μp' .+ σp'.*(randn(l,55).*(1 .- R) + x0*R)
	q = μq' .+ σq'.*(randn(l,55).*(1 .- R) + x0*R)
	
	v = zeros(0,n)
	for i in 1:l
		t = times[i]
#		vi = set_pq_and_solve(p[i,:],q[i,:],network_data)
		vi = set_pq_and_solve(Dict{String,Float64}("$k" => p[i,k] for k in 1:55), 
				      Dict{String,Float64}("$k" => q[i,k] for k in 1:55),
				      network_data
				      )
		v = [v; vi']
	end

	return p,q,v
end

p,q,v = gen_2p_correlated_time_series(30,32,.5)





