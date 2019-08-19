using LinearAlgebra

function objective1(Xs::Array{Float64,2}, Lm::Array{Float64,2}, dm::Array{Float64,1}, a::Array{Float64,1}, f::Array{Float64,1}, phi::Array{Float64,1}, dt::Float64, l1::Float64, l2::Float64)
	nn,T = size(Xs)
	n = Int(nn/2)

	err = Xs[n+1:2*n,2:end] - Xs[n+1:2*n,1:end-1] - dt * (-Lm*Xs[1:n,1:end-1] - diagm(0 => dm)*Xs[n+1:2*n,1:end-1] + repeat(a,1,T-1).*cos.(2*pi*dt*f*Array(1:T-1)' .+ repeat(phi,1,T-1)))
	
	return sum(err.^2) - l1 * sum(Lm .* (1 .- diagm(0 => ones(n)))) + l2 * sum(a)
end


function objective2(Xs::Array{Float64,2}, Lm::Array{Float64,2}, dm::Array{Float64,1}, a::Array{Float64,1}, f::Array{Float64,1}, phi::Array{Float64,1}, dt::Float64, l::Float64)
	nn,T = size(Xs)
	n = Int(nn/2)

	err = Xs[n+1:2*n,2*end] - Xs[n+1:2*n,1:end-1] - dt * (-Lm*Xs[1:n,1:end-1] - diagm(0 => dm)*Xs[n+1:2*n,1:end-1] + repeat(a,1,T-1).*cos.(2*pi*dt*f*Array(1:T-1)' .+ repeat(phi,1,T-1)))

	return sum(err.^2) + l*sum(a)
end


function objective3(Xs::Array{Float64,2}, Lm::Array{Float64,2}, dm::Array{Float64,1}, a::Array{Float64,1}, f::Array{Float64,1}, phi::Array{Float64,1}, dt::Float64, l::Float64)
	nn,T = size(Xs)
	n = Int(nn/2)

	err = zeros(n,T-1)

	for i in 1:n
		for t in 1:T-1
			err[i,t] = Xs[n+i,t+1] - Xs[n+i,t] - dt * (sum(-Lm[i,k]*Xs[k,t] for k = 1:n) - dm[i]*Xs[n+i,t] + a[i]*cos(2*pi*f[i]*t*dt + phi[i]))
		end
	end

	return sum(err[i,t]^2 for i = 1:n for t = 1:T-1)
end





