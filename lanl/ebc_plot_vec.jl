using PyPlot, LinearAlgebra

function ebc_plot_vec(v::Vector{Float64}, cols::String="plasma")
	ebc_plot_vec(v,vec(1:length(v)),cols)
end

function ebc_plot_vec(v::Vector{Float64}, ids::Vector{Int64}, cols::String="plasma")
	ixy = readdlm("data_ebc/coord.dat")
	idx = Int64.(ixy[:,1]) .+ 1
	i2i = Dict{Int64,Int64}()
	for i in 1:length(idx)
		i2i[idx[i]] = i
	end
	x = ixy[:,2]
	y = ixy[:,3]

	vmax = maximum(v)
	vmin = minimum(v)

	cmap = get_cmap(cols)

	for i in 1:length(ids)
		ii = ids[i]
		PyPlot.plot(x[i2i[ii]],y[i2i[ii]],"o",color=cmap((v[i]-vmin)/(vmax-vmin)))
	end

	return nothing
end








