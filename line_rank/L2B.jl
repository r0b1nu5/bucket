using SparseArrays

function L2B(L::Array{Float64,2})
	n = size(L)[1]
	
	B = Array{Float64,2}(undef,n,0)
	w = Array{Float64,1}(undef,0)
	for i in 1:n-1
		for j in i+1:n
			if L[i,j] != 0.0
				ed = zeros(n)
				ed[i] = 1.0
				ed[j] = -1.0
				B = [B ed]
				push!(w,-L[i,j])
			end
		end
	end
	
	return B,w
end

function L2B(L::SparseMatrixCSC{Float64,Int})
	n = size(L)[1]
	
	B = Array{Float64,2}(undef,n,0)
	w = Array{Float64,1}(undef,0)
	for i in 1:n-1
		for j in i+1:n
			if L[i,j] != 0.0
				ed = zeros(n)
				ed[i] = 1.0
				ed[j] = -1.0
				B = [B ed]
				push!(w,-L[i,j])
			end
		end
	end
	
	return sparse(B),w
end
