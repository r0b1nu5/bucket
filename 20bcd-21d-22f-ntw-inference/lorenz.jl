function lorenz(u0::Vector{Float64}, σ::Float64=10., ρ::Float64=28., β::Float64=8/3, h::Float64=.01, T::Int64=10000)
	u = copy(u0)
	us = copy(u0)
	dus = zeros(3,0)

	for t in 1:T
		k1 = f_lorenz(u,σ,ρ,β)
		k2 = f_lorenz(u+h*k1/2,σ,ρ,β)
		k3 = f_lorenz(u+h*k2/2,σ,ρ,β)
		k4 = f_lorenz(u+h*k3,σ,ρ,β)

		du = (k1 + 2*k2 + 2*k3 + k4)/6
		
		dus = [dus du]
		us = [us u+h*du]
		u = us[:,end]
	end

	return us,dus
end

function f_lorenz(u::Vector{Float64}, σ::Float64=10., ρ::Float64=28., β::Float64=8/3)
	x,y,z = u
	fx = σ*(y - x)		# -σx + σy
	fy = x*(ρ - z) - y 	# ρx - y - xz
	fz = x*y - β*z		# -βz + xy

	return [fx,fy,fz]
end

function f_lorenz(U::Matrix{Float64}, σ::Float64=10., ρ::Float64=28., β::Float64=8/3)
	n,T = size(U)
	fU = zeros(3,0)
	
	for t in 1:T
		fU = [fU f_lorenz(U[:,t],σ,ρ,β)]
	end

	return fU
end


