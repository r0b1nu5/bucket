# Compute the range of the kernel of Bout.

function kerrange(Bout::Array{Float64,2})
	n,m = size(Bout)

	d = Int.(vec(sum(Bout .!= 0.,dims=2)))

	A = zeros(m,0)

	for i in 1:length(d)
		ids = setdiff((Bout[i,:] .!= 0).*(1:m),[0,])

		a = zeros(m,d[i]-1)

		for k in 1:d[i]-1
			v1 = [ones(k);-k]
			v2 = v1/norm(v1)
			a[ids[1:k+1],k] = v2
		end

		A = [A a]
	end

	return A
end




