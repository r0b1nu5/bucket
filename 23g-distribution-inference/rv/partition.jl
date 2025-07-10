using Base.Iterators

function get_partition(list_meter::Vector{String}, list_pe::Vector{Vector{String}}, list_fs::Vector{Vector{String}})
	partition = Vector{Vector{String}}()
	remaining = copy(list_meter)
	p_done = String[]
	multi_p = false

	list_e,list_p = get_pe(list_pe)
	list_ef = setdiff(list_e,list_p) # List of children that are not parents

	for e in list_ef
		p = [pe[1] for pe in list_pe][[pe[2] for pe in list_pe] .== e]
		multi_p = (multi_p || length(p) > 1) # If multiple parents exist, we raise a flag and only keep the first one.

		if p[1] in p_done
			for part in partition
				if p[1] in part
					push!(part,e)
					remaining = setdiff(remaining,e)
				end
			end
		else
			push!(partition,[p[1],e])
			remaining = setdiff(remaining,[p[1],e])
		end
	end

#	list_fs_all = sort(union(collect(flatten(list_fs))))

	for fs in list_fs
		(fs[1] in list_p || fs[2] in list_p) && continue # If one of the siblings is a parent, we don't treat the pair.
		(fs[1] in list_e && fs[2] in list_e) && continue # If both of the siblings are children, they already have a partition.
		if !(fs[1] in remaining)
			for part in partition
				if fs[1] in part
					push!(part,fs[2])
					remaining = setdiff(remaining,fs[2])
				end
			end
		elseif !(fs[2] in remaining)
			for part in partition
				if fs[2] in part
					push!(part,fs[1])
					remaining = setdiff(remaining,fs[1])
				end
			end
		else
			push!(partition,["h";fs])
			remaining = setdiff(remaining,fs)
		end
	end

	for r in remaining
		push!(partition,[r,])
	end

	if multi_p
		@info "WARNING: Some children had multiple parents!"
	end

	return partition
end

function get_pe(list_pe::Vector{Vector{String}})
	list_e = Vector{String}()
	list_p = Vector{String}()

	for pe in list_pe
		push!(list_p,pe[1])
		push!(list_e,pe[2])
	end

	return list_e,list_p
end



