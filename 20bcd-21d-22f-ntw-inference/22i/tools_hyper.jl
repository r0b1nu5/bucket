function get_adj_3rd(A2,A3)
	n = size(A2)[1]

	adj = zeros(n,Int64(n*(n-1)/2))
	table = Dict{Int64,Tuple{Int64,Int64}}()

	for i in 1:n
		c = 0
		for n1 in 1:n-1
			for n2 in n1+1:n
				c += 1
				if i == n1
					adj[i,c] = A2[i,n2]
				elseif i == n2
					adj[i,c] = A2[i,n1]
				else
					adj[i,c] = A3[i,n1,n2]
				end
				table[c] = (n1,n2)
			end
		end
	end

	return adj,table
end

function adj2As(adj::Union{BitMatrix,Matrix{Float64}})
	n = size(adj)[1]

	A2 = zeros(n,n)
	A3 = zeros(n,n,n)

	for i in 1:n
		c = 0
		for n1 in 1:n-1
			for n2 in n1+1:n
				c += 1
				if i == n1
					A2[i,n2] = adj[i,c]
				elseif i == n2
					A2[i,n1] = adj[i,c]
				else
					A3[i,n1,n2] = adj[i,c]
				end
			end
		end
	end
	
	return A2,A3
end





