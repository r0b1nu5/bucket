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

function rm_line_inc(B::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int}},w::Array{Float64,1},l::Tuple{Int,Int})
	n,m = size(B)
	
	x,y,z = findnz(B)
	idx = setdiff(((x[2*(1:m).-1] .== l[1]).*(x[2*(1:m)] .== l[2]) + (x[2*(1:m).-1] .== l[2]).*(x[2*(1:m)] .== l[1])).*(1:m),[0,])
	
	if length(idx) == 0
		@info "The line does not exist..."
	elseif length(idx) > 1
		@info "Multiple lines ???"
		Br = B[:,setdiff(1:m,idx)]
		wr = w[setdiff(1:m,idx)]
	else
		Br = B[:,setdiff(1:m,idx)]
		wr = w[setdiff(1:m,idx)]
	end
		
	return Br,wr
end






