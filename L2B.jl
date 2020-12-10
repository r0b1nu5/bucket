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

function L2B(L::SparseMatrixCSC{Float64,Int},tol::Float64=1e-8)
	n = size(L)[1]
	
	I,J,V = findnz(L)
	mm = length(I)
	
	m = 0 
	IB = Array{Int64,1}()
	JB = Array{Int64,1}()
	VB = Array{Flaot64,1}()
	w = Array{Float64,1}()

	for k in 1:mm
		i = I[k]
		j = J[k]
		v = V[k]
		if i < j && abs(v) > tol
			m += 1
			push!(IB,i)
			push!(JB,m)
			push!(VB,1.)
			push!(IB,j)
			push!(JB,m)
			push!(VB,-1.)
			push!(w,v)
		end
	end

	B = sparse(IB,JB,VB)
	Bt = sparse(JB,IB,VB)

	return B,w,Bt
end

function L2B_ij(L::Array{Float64,2}, i0::Int64, j0::Int64)
	n = size(L)[1]

	if i0 > j0
		k = i0
		i0 = copy(j0)
		j0 = copy(i0)
	end

	B = Array{Float64,2}(undef,n,0)
	w = Array{Float64,1}(undef,0)
	c = 0
	eij = 0
	for i in 1:n-1
		for j in i+1:n
			if L[i,j] != 0.0 && i == i0 && j == j0
				c += 1
				eij = copy(c)
				ed = zeros(n)
				ed[i] = 1.0
				ed[j] = -1.0
				B = [B ed]
				push!(w,-L[i,j])
			elseif L[i,j] != 0.0
				c += 1
				ed = zeros(n)
				ed[i] = 1.0
				ed[j] = -1.0
				B = [B ed]
				push!(w,-L[i,j])
			end
		end
	end

	return B,w,eij
end


				 
