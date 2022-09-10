using LinearAlgebra 

function kf(L::Array{Float64,2},m::Int=1)
	@info("Computing Kfm...")
	
	n = size(L)[1]
	ls = sort(eigvals(L))[2:end]
	
	return n*sum(ls.^(-m))
end



