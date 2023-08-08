using PyPlot

function plot_A(A::Matrix{Float64},x::Vector{Float64},y::Vector{Float64})
	n = length(x)
	for i in 1:n
		for j in i+1:n
			if A[i,j] > 1e-10
				PyPlot.plot(x[[i,j]],y[[i,j]],"-k")
			end
		end
		PyPlot.plot(x[i],y[i],"ok")
	end
end

# Vertex weights are in v and are normalize
function plot_vscale(A::Matrix{Float64},x::Vector{Float64},y::Vector{Float64},v::Vector{Float64},cm::String="plasma";cb::Bool=false,cbl::String="")
	n = length(x)
	
	cmap = get_cmap(cm)

	for i in 1:n
		for j in i+1:n
			if A[i,j] > 1e-10
				PyPlot.plot(x[[i,j]],y[[i,j]],"-k")
			end
		end
		PyPlot.plot(x[i],y[i],"o",color=cmap(v[i]))
	end
end

# Edge weights are in A and are normalized.
function plot_escale(A::Matrix{Float64},x::Vector{Float64},y::Vector{Float64},cm::String="plasma";cb::Bool=false,cbl::String="")
	n = length(x)

	cmap = get_cmap(cm)

	for i in 1:n
		for j in i+1:n
			if abs(A[i,j]) > 1e-10
				PyPlot.plot(x[[i,j]],y[[i,j]],color=cmap(A[i,j]))
			end
		end
		PyPlot.plot(x[i],y[i],"ok")
	end
end

function get_list(v::Vector{Float64},names::Vector{String},n::Int64)
	V = sortslices([v names],dims=1,rev=true)
	V1 = Float64.(V[:,1])
	V2 = String.(V[:,2])
	str = ""
	for i in 1:n
		str *= V2[i]*": $(V1[i])"*'\n'
	end
	return str
end

function normalized(v::Union{Vector{Float64},Vector{Int64}})
	return (v .- minimum(v))./(maximum(v) - minimum(v))
end



