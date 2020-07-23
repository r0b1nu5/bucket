using LinearAlgebra, LightGraphs


function small_world(n::Int64, k::Int64, p::Float64)
	L = diagm(0 => 2*k*ones(n))

	for i in 1:k
		L += diagm(i => -ones(n-i)) + diagm(-i => -ones(n-i)) + diagm(n-i => -ones(i)) + diagm(i-n => -ones(i))
	end

	for i in 1:k
		for j in 1:n
			if rand() < p
				l = rand(Int.(setdiff((1:n).*(L[j,:] .== 0.),[0,])))
				L[[j,mod(j+k-1,n)+1],[j,mod(j+k-1,n)+1]] -= [1. -1.;-1. 1.]
				L[[j,l],[j,l]] += [1. -1;-1. 1.]
			end
		end
	end

	return L
end


function clustering(L::Array{Float64,2})
	n = size(L)[1]

	trip = 0
	ctrip = 0
	
	for i in 1:n-2
		for j in i+1:n-1
			for k in j+1:n
				if L[i,j] < 0 && L[i,k] < 0
					trip += 1
					if L[j,k] < 0
						ctrip += 1
					end
				end
			end
		end
	end

	return ctrip/trip
end


function avg_path_l(L::Array{Float64,2})
	n = size(L)[1]
	
	g = SimpleGraph()

	add_vertices!(g,n)

	for i in 1:n-1
		for j in i+1:n
			if L[i,j] < 0.
				add_edge!(g,i,j)
			end
		end
	end

	Ds = zeros(n,n)
	for i in 1:n
		Ds[:,i] = dijkstra_shortest_paths(g,i).dists
	end

	a = sum(Ds)*2/(n*(n-1))

	return a
end




