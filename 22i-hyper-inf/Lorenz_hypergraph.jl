
σ = 10.0
ρ = 28.0
β = 8.0 / 3.0
k2 = 1.0 # Pairwise coupling strength
k3 = 1.0 # Three-body coupling strength

function lorenz(t, xyz, σ, ρ, β)
	x,y,z = xyz
	dx = σ*(y - x)
	dy = x*(ρ - z) - y
	dz = x*y - β*z
	return [dx,dy,dz]
end

function coupled_lorenz(t, xyz::Vector{Float64}, σ, ρ, β, links, triangles, k2,k3)
	N = Int64(length(xyz)/3)
	dx = zeros(length(xyz))

	for i in 1:N
		xi,yi,zi = xyz[(i-1)*3 .+ (1:3)]
		dx[(i-1)*3 .+ (1:3)] += lorenz(t,[xi,yi,zi],σ,ρ,β)
	end

	for ij in links
		i,j = ij
		xi,yi,zi = xyz[(i-1)*3 .+ (1:3)]
		xj,yj,zj = xyz[(j-1)*3 .+ (1:3)]
		dx[(i-1)*3+1] += k2*(xj-xi)
		dx[(j-1)*3+1] += k2*(xi-xj)
	end

	for ijk in triangles
		i,j,k = ijk
		xi,yi,zi = xyz[(i-1)*3 .+ (1:3)]
		xj,yj,zj = xyz[(j-1)*3 .+ (1:3)]
		xk,yk,zk = xyz[(k-1)*3 .+ (1:3)]
		dx[(i-1)*3+1] += k3*(xk*xj^2-xi^2)
		dx[(j-1)*3+1] += k3*(xk*xi^2-xj^2)
		dx[(k-1)*3+1] += k3*(xj*xi^2-xk^2)
	end

	return dx
end

function coupled_lorenz(t, xyz::Matrix{Float64}, σ, ρ, β, links, triangles, k2,k3)
	N,T = size(xyz)
	y = zeros(N,0)
	for i in 1:T
		y = [y coupled_lorenz(t,xyz[:,i],σ,ρ,β,links,triangles,k2,k3)]
	end

	return y
end





