using Statistics, LinearAlgebra, PyPlot

# X: Data matrix. Each row is one time step, each column is one type of data.
# c0: Initial centers of the groups. Each row is one center
function my_kmeans(X::Matrix{Float64}, c0::Matrix{Float64}, γ::Vector{Float64}, max_iter::Int64=1000)
	if size(X)[2] != size(c0)[2]
		@info "Sizes don't match..."
	end

	k,m = size(c0)
	n,m = size(X)

	c = c0
	g1 = zeros(Int64,n)
	g2 = ones(Int64,n)

	iter = 0
	while g1 != g2 && iter < max_iter
		iter += 1

		D = dist2c(X,c,γ)
		g1 = g2
		g2 = [findmin(D[i,:])[2] for i in 1:n]
		for i in 1:k
			c[i,:] = vec(mean(X[g2 .== i,:],dims=1))
		end
	end

	return g2, c
end

function my_kmeans(X::Matrix{Float64}, c0::Matrix{Float64}, max_iter::Int64=1000)
	return my_kmeans(X,c0,ones(size(c0)[2]),max_iter)
end


function dist2c(X::Matrix{Float64}, c::Matrix{Float64}, γ::Vector{Float64})
	if size(X)[2] != size(c0)[2] || size(X)[2] != length(γ) || size(c)[2] != length(γ)
		@info "Sizes don't match..."
	end

	n = size(X)[1]
	k = size(c)[1]

	D = zeros(n,k)

	for i in 1:n
		for j in 1:k
			D[i,j] = norm(γ.*(X[i,:]-c[j,:]))
		end
	end

	return D
end

function dist2c(X::Matrix{Float64}, c::Matrix{Float64})
	return dist2c(X,c,ones(size(c)[2]))
end


function plot_grps(X::Matrix{Float64}, g::Vector{Int64}, dims::Vector{Int64})
	if maximum(dims) > size(X)[2]
		@info "Not enough dimensions to match dims. Aborting."
		return nothing
	end
	if minimum(dims) < 1
		@info "Ill-defined index in dims. Aborting."
		return nothing
	end

	if length(dims) < 1
		@info "Empty plot. Aborting."
		return nothing
	elseif length(dims) == 1
		c = [mean(X[g .== i,dims]) for i in 1:maximum(g)]
		for i in 1:maximum(g)
			PyPlot.plot(X[g .== i,dims[1]],i*ones(length(g[g .== i])),".",color="C$(i-1)")
			PyPlot.plot([c[i],c[i]],[-1,1] .+ i,"--",color="C$(i-1)")
		end
	elseif length(dims) == 2
		c = [vec(mean(X[g .== i,dims],dims=1)) for i in 1:maximum(g)]
		for i in 1:maximum(g)
			PyPlot.plot(X[g .== i,dims[1]],X[g .== i,dims[2]],".",color="C$(i-1)")
			PyPlot.plot(c[i][1],c[i][2],"s",color="C$(i-1)")
		end
	elseif length(dims) == 3
		c = [vec(mean(X[g .== i,dims],dims=1)) for i in 1:maximum(g)]
		for i in 1:maximum(g)
			PyPlot.plot3D(X[g .== i,dims[1]],X[g .== i,dims[2]],X[g .== i,dims[3]],".",color="C$(i-1)")
			PyPlot.plot3D(c[i][1],c[i][2],c[i][3],"s",color="C$(i-1)")
		end
	else
		@info "Too many dimensions to show. Keeping only the first two."
		plot_grps(X,g,dims[1:2])
	end

	return nothing
end

function plot_grps(X::Matrix{Float64}, g::Vector{Int64})
	plot_grps(X,g,[1,2])
end







