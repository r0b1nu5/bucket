using LinearAlgebra

function objective1(Xs::Array{Float64,2}, Lm::Array{Float64,2}, dm::Array{Float64,1}, a::Array{Float64,1}, f::Array{Float64,1}, phi::Array{Float64,1}, dt::Float64, l::Float64)
	nn,T = size(Xs)
	n = Int(nn/2)

	err = Xs[n+1:2*n,2:end] - Xs[n+1:2*n,1:end-1] - dt * (Lm*Xs[1:n,1:end-1] - diagm(0 => dm)*Xs[n+1:2*n,1:end-1] - repeat(a,1,T-1).*cos.(2*pi*dt*f*Array(1:T-1)' .+ repeat(phi,1,T-1)))
	
	return sum(err.^2) - l*sum(Lm .* (1 .- diagm(0 => ones(n))))
end


function objective2(Xs::Array{Float64,2}, Lm::Array{Float64,2}, dm::Array{Float64,1}, a::Array{Float64,1}, f::Array{Float64,1}, phi::Array{Float64,1}, dt::Float64, l::Float64)
	nn,T = size(Xs)
	n = Int(nn/2)

	err = Xs[n+1:2*n,2*end] - Xs[n+1:2*n,1:end-1] - dt * (Lm*Xs[1:n,1:end-1] - diagm(0 => dm)*Xs[n+1:2*n,1:end-1] - repeat(a,1,T-1).*cos.(2*pi*dt*f*Array(1:T-1)' .+ repeat(phi,1,T-1)))

	return sum(err.^2) + l*sum(a)
end




