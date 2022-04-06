
function rand_idx(n::Int64, k::Int64)
	ids = Vector(1:n)
	idx = Vector{Int64}()

	for i in 1:k
		id = rand(ids)
		push!(idx,id)
		ids = setdiff(ids,[id,])
	end

	return idx
end




