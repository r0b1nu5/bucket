using LinearAlgebra, DelimitedFiles, Distributions

function generate_time_series(L::Array{Float64,2}, d::Array{Float64,1}, m::Array{Float64,1}, T::Int64, dt::Float64)
	n = size(L)[1]
	
	Mi = diagm(0 => 1 ./ m)
	D = diagm(0 => d)
	
	X = zeros(2*n)
#	X = [2*pi*rand(n,1);zeros(n,1)]
	
	t = 0
	
	while t < T
		t += 1
		
		xi = rand(Normal(0,1),n)
		
		X = [X (X[:,end] + dt * [X[n+1:2*n,end] ; -Mi*D*X[n+1:2*n,end] - Mi*L*X[1:n,end] + Mi*xi./sqrt(dt)])]
	end
	
	return X
end








