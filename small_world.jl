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


function resistive_world(N::Int64, k::Int64)
	d = round(Int,k/2)

	n = k+1

	L = 2*diagm(0 => ones(n)) - diagm(1 => ones(n-1)) - diagm(-1 => ones(n-1)) - diagm(n-1 => ones(1)) - diagm(1-n => ones(1))
	Ld = pinv(L)

	while n < N
		n += 1
		if n%10 == 0
			@info "n = $n"
		end

		i = rand(1:n-1)

		L = [L zeros(n-1);zeros(1,n)]
		Ld = [Ld zeros(n-1);zeros(1,n)]
		
		ei = zeros(n)
		ei[[i,n]] = [1.,-1.]
		
		L += ei*ei'
		Ld = Ld - (Ld*ei*ei'*Ld)/(1 + ei'*Ld*ei)

		Omn = Ld[n,n] .+ diag(Ld) - 2*Ld[n,:]

		tolerate = setdiff(Array(1:n-1),[i,])
		for l in 2:k
			p = 1 ./Omn[tolerate]
			p ./= sum(p)
			c = [sum(p[1:j]) for j in 0:length(tolerate)]

			x = rand()
			new_id_t = setdiff((x .> c).*(1:length(tolerate)+1),[0,])[end]
			new_id = tolerate[new_id_t]

			tolerate = setdiff(tolerate,[new_id,])

			ei = zeros(n)
			ei[[n,new_id]] = [1.,-1.]

			L += ei*ei'
			Ld = Ld - (Ld*ei*ei'*Ld)/(1 + ei'*Ld*ei)
		end

	end

	return L
end




