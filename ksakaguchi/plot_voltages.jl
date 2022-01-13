using PyPlot

function plot_voltages(v::Vector{ComplexF64}, ids::Vector{Int64}, col::Any, vmax::Float64=1.1, vmin::Float64=.9)
	if vmax > 0.
		PyPlot.plot(vmax*cos.(LinRange(0,2π,200)),vmax*sin.(LinRange(0,2π,200)),"--",color="C7",linewidth=.5)
		PyPlot.plot(vmin*cos.(LinRange(0,2π,200)),vmin*sin.(LinRange(0,2π,200)),"--",color="C7",linewidth=.5)
	end
	
	if length(ids) == 0
		PyPlot.plot(real.(v),imag.(v),"o",color=col)
	else
		PyPlot.plot(real.(v[ids]),imag.(v[ids]),"-o",color=col)
		for i in ids
			PyPlot.text(real(v[i]),imag(v[i]),"$i")
		end
	end

	return nothing
end

function plot_voltages(v::Vector{ComplexF64}, ids::Vector{Vector{Int64}}, col::Vector{Any}, vmax::Float64=1.1, vmin::Float64=.9)
	if vmax > 0.
		PyPlot.plot(vmax*cos.(LinRange(0,2π,200)),vmax*sin.(LinRange(0,2π,200)),"--",color="C7",linewidth=.5)
		PyPlot.plot(vmin*cos.(LinRange(0,2π,200)),vmin*sin.(LinRange(0,2π,200)),"--",color="C7",linewidth=.5)
	end

	c = 0
	for id in ids
		c += 1
		plot_voltages(v,id,col[c]-1.)
	end

	return nothing
end

function plot_voltages(vs::Vector{Vector{ComplexF64}}, ids::Vector{Int64}, col::Vector{Any}, vmax::Float64=1.1, vmin::Float64=.9)
	if vmax > 0.
		PyPlot.plot(vmax*cos.(LinRange(0,2π,200)),vmax*sin.(LinRange(0,2π,200)),"--",color="C7",linewidth=.5)
		PyPlot.plot(vmin*cos.(LinRange(0,2π,200)),vmin*sin.(LinRange(0,2π,200)),"--",color="C7",linewidth=.5)
	end

	c = 0
	for v in vs
		c += 1
		plot_voltages(v,ids,col[c],-1.)
	end

	return nothing
end




