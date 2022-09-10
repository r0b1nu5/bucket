function smooth85(xs::Vector{Float64}, id0::Int64=8, did::Int64=10)
	ids = vec(id0:did:length(xs))

	corr = (xs[ids[1:end-1].-1] + xs[ids[1:end-1].+1])./2

	ys = xs
	ys[ids[1:end-1]] = corr

	return ys
end




