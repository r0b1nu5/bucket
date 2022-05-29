using PyPlot

function plot_neighbors(id::Int64, L::Matrix{Float64}, X::Array{Float64,3}, α::Float64=.1)
	l1,l2,l3 = size(X)
	x = zeros(0,l3)
	for i in 1:l1
		x = [x;X[i,:,:]]
	end

	plot_neighbors(id,L,x,α)
end

function plot_neighbors(id::Int64, L::Matrix{Float64}, xs::Matrix{Float64}, X::Array{Float64,3}, α::Float64=.1)
	l1,l2,l3 = size(X)
	x = zeros(0,l3)
	for i in 1:l1
		x = [x;X[i,:,:]]
	end

	plot_neighbors(id,L,xs,x,α)
end

function plot_neighbors(id::Int64, L::Matrix{Float64}, X::Matrix{Float64}, α::Float64=.1)
	k,l = size(X)
	n = size(L)[1]
	
	ids = setdiff((abs.(L[id,:]) .> 1e-8).*(1:n),[0,])
	
	figure("traj",(10.,10.))
	PyPlot.plot3D(X[:,1],X[:,2],X[:,3],"o",color="C0",alpha=α)
	PyPlot.plot3D(X[ids,1],X[ids,2],X[ids,3],"o",color="C1")
end

function plot_neighbors(id::Int64, L::Matrix{Float64}, xs::Matrix{Float64}, X::Matrix{Float64}, α::Float64=.1)
	k,l = size(X)
	n = size(L)[1]
	
	ids = setdiff((abs.(L[id,:]) .> 1e-8).*(1:n),[0,])

	figure("traj",(10.,10.))
	PyPlot.plot3D(X[:,1],X[:,2],X[:,3],"o",color="C0",alpha=α)
	PyPlot.plot3D(X[ids,1],X[ids,2],X[ids,3],"o",color="C1")
	PyPlot.plot(xs[1,:],xs[2,:],xs[3,:],color="C3")
end


function gen_grid(res::Int64)
	x = Vector{Float64}()
	y = Vector{Float64}()
	z = Vector{Float64}()

	δs = LinRange(0.,1.,res)

	for i in 1:res
		for j in 1:res-i+1
			push!(x,δs[i])
			push!(y,δs[j])
			push!(z,1-δs[i]-δs[j])
		end
	end

	return x,y,z
end

function plot_radius(x0::Vector{Float64}, ϵ::Float64, res::Int64, col::String="C2", α::Float64=1.)
	x,y,z = gen_grid(res)

	rx = Vector{Float64}()
	ry = Vector{Float64}()
	rz = Vector{Float64}()
	
	for i in 1:length(x)
		if norm(x0 - [x[i],y[i],z[i]],1) < ϵ
			push!(rx,x[i])
			push!(ry,y[i])
			push!(rz,z[i])
		end
	end
	
	PyPlot.plot3D(rx,ry,rz,".",color=col,alpha=α)
end


