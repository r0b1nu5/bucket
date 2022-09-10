using LinearAlgebra

# Creates an oriented cactus graph

function cactus(cyc::Array{Tuple{Int64,Int64}})
	ns = Array{Int64,1}()
	js = Array{Int64,1}()

	for c in cyc
		push!(ns,c[1])
		push!(js,c[2])
	end

	N = sum(ns) - length(ns) + 1

	L = zeros(N,N)
	L[1:ns[1],1:ns[1]] += diagm(0 => ones(ns[1])) - diagm(1 => ones(ns[1]-1)) - diagm((1-ns[1]) => ones(1))
	nt = ns[1]

	for i in 2:length(ns)
		n = ns[i]
		j = js[i]

		ids = [j;(nt+1):(nt+n-1)]

		L[ids,ids] += diagm(0 => ones(n)) - diagm(1 => ones(n-1)) - diagm((1-n) => ones(1))

		nt += (n-1)
	end

	return L
end
	



