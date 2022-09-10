function histo(list::Array{Float64,1},bins::Array{Float64,1})
	h = Array{Float64,1}()
	for i in 1:(length(bins)-1)
		push!(h,sum(bins[i] .<= list .< bins[i+1]))
	end
	return h
end

function histo(list::Array{Int64,1},bins::Array{Float64,1})
	h = Array{Float64,1}()
	for i in 1:(length(bins)-1)
		push!(h,sum(bins[i] .<= list .< bins[i+1]))
	end
	return h
end


