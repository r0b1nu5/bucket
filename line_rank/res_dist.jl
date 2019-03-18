using LinearAlgebra,SparseArrays

function res_dist(L::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int64}})
	n = size(L)[1]
	
	Gai = pinv(Array(L))
	
	Om = zeros(n,n)
	for i in 1:n-1
		for j in i+1:n
			Om[i,j] = Gai[i,i] + Gai[j,j] - 2*Gai[i,j]
			Om[j,i] = Om[i,j]
		end
	end
	
	return Om
end



