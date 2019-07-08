using LinearAlgebra, DelimitedFiles, Distributions

function generate_time_series(ntw::String, L::Array{Float64,2}, m::Array{Float64,1}, d::Array{Float64,1}, T::Int64, dt::Float64, sig::Array{Float64,1})
	n = size(L)[1]
	
	Mi = diagm(0 => 1 ./ m)
	D = diagm(0 => d)
	
	Ad = [zeros(n,n) diagm(0 => ones(n)); -Mi*L -Mi*D]
	A = Ad*dt + diagm(0 => ones(2*n))
	
	B = diagm(0 => sqrt(dt)*sig./m)
	
	X = zeros(2*n)
#	X = [2*pi*rand(n,1);zeros(n,1)]
	
	t = 0
	
	while t < T
		t += 1
		
		xi = rand(Normal(0,1),n)
		
		X = [X (A*X[:,end] + [zeros(n);B*xi])]
	end
	
	writedlm("data/"*ntw*"_$(T)_$(dt).csv",X,',')
	
	return X
end








