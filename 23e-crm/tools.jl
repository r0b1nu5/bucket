using PyPlot

function plot_A(A::Matrix{Float64},x::Vector{Float64},y::Vector{Float64})
	n = length(x)
	for i in 1:n
		for j in i+1:n
			if A[i,j] > 1e-10
				PyPlot.plot(x[[i,j]],y[[i,j]],"-k",zorder=-1)
			end
		end
	end
	PyPlot.scatter(x,y;s=70*ones(length(x)),color="k",zorder=1)

	dx = maximum(x) - minimum(x)
	xmin = minimum(x) - .1*dx
	xmax = maximum(x) + .1*dx
	dy = maximum(y) - minimum(y)
	ymin = minimum(y) - .1*dy
	ymax = maximum(y) + .1*dy

	axis([xmin,xmax,ymin,ymax])
end

# Vertex weights are in v
function plot_vscale(A::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64}, v::Union{Vector{Int64},Vector{Float64}}, cm::String="plasma"; cb::Bool=false, cbl::String="")
	n = length(x)
	
	cmap = get_cmap(cm)

	for i in 1:n
		for j in i+1:n
			if A[i,j] > 1e-10
				PyPlot.plot(x[[i,j]],y[[i,j]],"-k",zorder=-1)
			end
		end
	end
	PyPlot.scatter(x,y;s=70*ones(length(x)),c=v,cmap=cm,zorder=1,edgecolors="k")

	if cb
		colorbar(label=cbl)
	end

	dx = maximum(x) - minimum(x)
	xmin = minimum(x) - .1*dx
	xmax = maximum(x) + .1*dx
	dy = maximum(y) - minimum(y)
	ymin = minimum(y) - .1*dy
	ymax = maximum(y) + .1*dy

	axis([xmin,xmax,ymin,ymax])
	xticks([])
	yticks([])
end

function plot_vscale(A::Matrix{Float64}, Add::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64}, v::Union{Vector{Int64},Vector{Float64}}, cm::String="plasma"; cb::Bool=false, cbl::String="")
	n = length(x)

	for i in 1:n-1
		for j in i+1:n
			if Add[i,j] > 1e-10
				PyPlot.plot(x[[i,j]],y[[i,j]],"--r",linewidth=2,zorder=-1)
			end
		end
	end

	plot_vscale(A,x,y,v,cm,cb=cb,cbl=cbl)
end

# Edge weights are in A
function plot_escale(A::Matrix{Float64},x::Vector{Float64},y::Vector{Float64},cm::String="plasma";cb::Bool=false,cbl::String="")
	n = length(x)

	cmap = get_cmap(cm)
	dA = maximum(A) - minimum(A)

	for i in 1:n
		for j in i+1:n
			if abs(A[i,j]) > 1e-10
				PyPlot.plot(x[[i,j]],y[[i,j]],color=cmap((A[i,j] - minimum(A))/dA),lw=2)
			end
		end
		PyPlot.plot(x[i],y[i],"ok")
	end

	if cb
		V = union(vec(A))
		PyPlot.scatter(-1000*ones(length(V)),-1000*ones(length(V));c=V,cmap=cm)
		colorbar(label=cbl)
	end
	
	dx = maximum(x) - minimum(x)
	xmin = minimum(x) - .1*dx
	xmax = maximum(x) + .1*dx
	dy = maximum(y) - minimum(y)
	ymin = minimum(y) - .1*dy
	ymax = maximum(y) + .1*dy

	axis([xmin,xmax,ymin,ymax])
	xticks([])
	yticks([])
end

function plot_vescale(A::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64}, v::Vector{Float64}, cmv::String="plasma", cme::String="plasma"; cbv::Bool=false, cbe::Bool=false, cbvl::String="", cbel::String="", fmax::Float64=-1.)
	n = length(x)

	fm = fmax
	if fmax < 0.
		fm = maximum(abs.(A))
	end

	cmapv = get_cmap(cmv)
	cmape = get_cmap(cme)

	for i in 1:n-1
		for j in i+1:n
			if abs(A[i,j]) > 1e-10
				PyPlot.plot(x[[i,j]],y[[i,j]],color=cmape(clamp(A[i,j]/fm,0.,1.)),zorder=-1,lw=2)
			end
		end
	end

	if cbe
		PyPlot.scatter([-1000,-1000],[-1000,-1000];c=[0.,fm],cmap=cme)
		#V = union(vec(A))
		#PyPlot.scatter(-1000*ones(length(V)),-1000*ones(length(V));c=V,cmap=cme)
		colorbar(label=cbel)
	end

#	PyPlot.scatter(x,y;s=70*ones(length(x)),c=v,cmap=cmv,zorder=1,edgecolors="k")
# #=
	X = [x;1000;1000]
	Y = [y;1000;1000]
	V = [v;maximum(abs.(v))*[-1,1]]
	PyPlot.scatter(X,Y;s=70*ones(length(X)),c=V,cmap=cmv,zorder=1,edgecolors="k")
# =#
 #=
	X = [x;1000;1000]
	Y = [y;1000;1000]
	V = [sign.(v).*log.(abs.(v));log.(maximum(abs.(v)))*[-1,1]]
	PyPlot.scatter(X,Y;s=70*ones(length(X)),c=V,cmap=cmv,zorder=1,edgecolors="k")
# =#
	if cbv
		colorbar(label=cbvl)
	end

	dx = maximum(x) - minimum(x)
	xmin = minimum(x) - .1*dx
	xmax = maximum(x) + .1*dx
	dy = maximum(y) - minimum(y)
	ymin = minimum(y) - .1*dy
	ymax = maximum(y) + .1*dy

	axis([xmin,xmax,ymin,ymax])
	xticks([])
	yticks([])
end

function plot_vescale(A::Matrix{Float64}, Add::Matrix{Float64}, x::Vector{Float64}, y::Vector{Float64}, v::Union{Vector{Int64},Vector{Float64}}, cmv::String="plasma", cme::String="plasma"; cbv::Bool=false, cbvl::String="", cbe::Bool=false, cbel::String="", fmax::Float64=-1.)
	n = length(x)
	fm = fmax
	if fmax < 0.
		maximum(abs.(A+Add))
	end

	cmape = get_cmap(cme)

	for i in 1:n-1
		for j in i+1:n
			if Add[i,j] > 1e-10
				PyPlot.plot(x[[i,j]],y[[i,j]],"--",color=cmape(clamp(Add[i,j]/fm,0.,1.)),linewidth=2,zorder=-1)
			end
		end
	end

	plot_vescale(A,x,y,v,cmv,cme,cbv=cbv,cbvl=cbvl,cbe=cbe,cbel=cbel,fmax=fm)
end

function get_list(v::Union{Vector{Int64},Vector{Float64}},names::Vector{String},n::Int64)
	V = sortslices([v names collect(1:length(v))],dims=1,rev=true)
	V1 = Float64.(V[:,1])
	V2 = String.(V[:,2])
	str = ""
	for i in 1:n
		str *= V2[i]*": $(V1[i])"*'\n'
	end
	return str,V2[1:n],V1[1:n],V[:,3]
end

function normalized(v::Union{Vector{Float64},Vector{Int64}})
	return (v .- minimum(v))./(maximum(v) - minimum(v))
end

function plot_ch(col::String="k", lw::Union{Float64,Int64}=1)
	xy = readdlm("xy-ch.csv",',')
	x = xy[:,1]
	y = xy[:,2]

	PyPlot.plot(x,y,color=col,linewidth=lw,zorder=-2)
end

function get_measure(m::String, g::SimpleGraph{Int64})
	if m == "bc"
		return betweenness_centrality(g)
	elseif m == "cc"
		return closeness_centrality(g)
	elseif m == "dc"
		return degree_centrality(g)
	elseif m == "de"
		return degree(g)
	elseif m == "ec"
		return eigenvector_centrality(g)
	elseif m == "pr"
		return pagerank(g)
	elseif m == "rc"
		return radiality_centrality(g)
	elseif m == "sc"
		return stress_centrality(g)
	else
		@info "Mesure inconnue (pour le moment)."
		return zeros(nv(g))
	end
end





