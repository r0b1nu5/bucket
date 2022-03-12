using PyPlot, LinearAlgebra, DelimitedFiles, SparseArrays

# !!! Directed graph !!!
function Bata(n::Int64)
	m = Int64(n*(n-1)/2)

	I = diagm(0 => ones(n-1))

	B = [ones(1,n-1);-I]
	for i in 2:n
		B = [B [-I[1:i-1,:];ones(1,n-1);-I[i:n-1,:]]]
	end

	return B
end

function spBata(n::Int64)
	I = Vector{Int64}()
	J = Vector{Int64}()
	V = Vector{Float64}()
	
	c = 0
	for i in 1:n
		for j in 1:n
			if i!=j
				c += 1
				push!(I,i)
				push!(J,c)
				push!(V,1.)
				push!(I,j)
				push!(J,c)
				push!(V,-1.)
			end
		end
	end
	
	return sparse(I,J,V)
end


# Computes the incidence matrix (B) and the edge weights (w) from the Laplacian matrix (l).
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
	VB = Array{Float64,1}()
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
			push!(w,-v)
		end
	end

	B = sparse(IB,JB,VB)
	Bt = sparse(JB,IB,VB)

	return B,w,Bt
end



