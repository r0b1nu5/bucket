using SparseArrays

function rm_line(L::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}},l::Tuple{Int,Int},tolerance::Float64=1e-8)
	nL = copy(L)
	if abs(L[l[1],l[2]]) < tolerance
		@info("Line $l does not exist")
	else
		nL[[l[1],l[2]],[l[1],l[2]]] -= L[l[1],l[2]]*[-1 1;1 -1]
	end
	
	return nL
end



