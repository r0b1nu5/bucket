using LightGraphs, SimpleWeightedGraphs, SparseArrays

function res_dist(L::Array{Float64,2})
	n = size(L)[1]
	Ld = pinv(L)

	return repeat(diag(Ld),1,n) + repeat(diag(Ld)',n,1) - Ld - Ld'
end

function res_dist(i::Int64, j::Int64, L::Array{Float64,2})
	Ω = res_dist(L)

	return Ω[i,j]
end

function geo_dist(L::Union{SparseMatrixCSC{Float64,Int64},Array{Float64,2}})
	n = size(L)[1]

	if issparse(L)
		A = spdiagm(0 => diag(L)) - L
	else
		A = diagm(0 => diag(L)) - L
	end

	g = SimpleWeightedGraph((A + A')/2)

	d = zeros(n,n)
	for i in 1:n
		d[:,i] = dijkstra_shortest_paths(g,i).dists
	end
	
	return d
end

function geo_dist(i::Int64, L::Union{SparseMatrixCSC{Float64,Int64},Array{Float64,2}})
	if issparse(L)
		A = spdiagm(0 => diag(L)) - L
	else
		A = diagm(0 => diag(L)) - L
	end

	g = SimpleWeightedGraph((A + A')/2)

	return dijkstra_shortest_paths(g,i).dists
end


function geo_dist(i::Int64, j::Int64, L::Union{SparseMatrixCSC{Float64,Int64},Array{Float64,2}})
	return geo_dist(i,L)[j]
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
