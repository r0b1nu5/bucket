using LinearAlgebra

#################################################################
# Get sensitivity matrices 
# Network matrices (G, B, Y) are assumed "Laplacian" (i.e., diagonal terms are minus the sum of the off-diagonal in the same row).
function get_Jpt(θ::Vector{Float64}, V::Vector{Float64}, B::Matrix{Float64}, G::Matrix{Float64})
	n = length(θ)
	Id = diagm(0 => ones(n))
	dI = 1 .- Id

	vs = V.*sin.(θ)
	vc = V.*cos.(θ)

	Jpt = dI.*(G.*(vs*vc' - vc*vs') - B.*(vc*vc' + vs*vs'))

	return Jpt - diagm(0 => vec(Jpt*ones(n)))
end

function get_Jpt(V::Vector{Complex{Float64}}, Y::Matrix{Complex{Float64}})
	return get_Jpt(angle.(V),abs.(V),imag.(Y),real.(Y))
end

function get_Jqt(θ::Vector{Float64}, V::Vector{Float64}, B::Matrix{Float64}, G::Matrix{Float64})
	n = length(θ)
	Id = diagm(0 => ones(n))
	dI = 1 .- Id

	vs = V.*sin.(θ)
	vc = V.*cos.(θ)

	Jqt = -dI.*(G.*(vc*vc' + vs*vs') + B.*(vs*vc' - vc*vs'))

	return Jqt - diagm(0 => vec(Jqt*ones(n)))
end

function get_Jqt(V::Vector{Complex{Float64}}, Y::Matrix{Complex{Float64}})
	return get_Jqt(angle.(V),abs.(V),imag.(Y),real.(Y))
end

function get_Jpv(θ::Vector{Float64}, V::Vector{Float64}, B::Matrix{Float64}, G::Matrix{Float64})
	n = length(θ)
	Id = diagm(0 => ones(n))
	dI = 1 .- Id

	s = sin.(θ)
	c = cos.(θ)
	vs = V.*s
	vc = V.*c

	Jpv = dI.*(G.*(vc*c' + vs*s') + B.*(vs*c' - vc*s'))

	return Jpv + diagm(0 => (2*V.*diag(G) + Jpv*ones(n)))
end

function get_Jpv(V::Vector{Complex{Float64}}, Y::Matrix{Complex{Float64}})
	return get_Jpv(angle.(V),abs.(V),imag.(Y),real.(Y))
end

function get_Jqv(θ::Vector{Float64}, V::Vector{Float64}, B::Matrix{Float64}, G::Matrix{Float64})
	n = length(θ)
	Id = diagm(0 => ones(n))
	dI = 1 .- Id

	s = sin.(θ)
	c = cos.(θ)
	vs = V.*s
	vc = V.*c

	Jqv = dI.*(G.*(vs*c' - vc*s') - B.*(vc*c' + vs*s'))

	return Jqv + diagm(0 => (-2*V.*diag(B) + Jqv*ones(n)))
end

function get_Jqv(V::Vector{Complex{Float64}}, Y::Matrix{Complex{Float64}})
	return get_Jqv(angle.(V),abs.(V),imag.(Y),real.(Y))
end

function get_J(θ::Vector{Float64}, V::Vector{Float64}, B::Matrix{Float64}, G::Matrix{Float64})
	n = length(θ)
	Id = diagm(0 => ones(n))
	dI = 1 .- Id

	s = sin.(θ)
	c = cos.(θ)
	vs = V.*s
	vc = V.*c

	Jpt = dI.*(G.*(vs*vc' - vc*vs') - B.*(vc*vc' + vs*vs'))
	Jpt -= diagm(0 => vec(Jpt*ones(n)))
	
	Jqt = -dI.*(G.*(vc*vc' + vs*vs') + B.*(vs*vc' - vc*vs'))
	Jqt -= diagm(0 => vec(Jqt*ones(n)))

	Jpv = dI.*(G.*(vc*c' + vs*s') + B.*(vs*c' - vc*s'))
	Jpv += diagm(0 => (2*V.*diag(G) + Jpv*ones(n)))
	
	Jqv = dI.*(G.*(vs*c' - vc*s') - B.*(vc*c' + vs*s'))
	Jqv += diagm(0 => (-2*V.*diag(B) + Jqv*ones(n)))

	return [Jpt Jpv;Jqt Jqv]
end

function get_J(V::Vector{Complex{Float64}}, Y::Matrix{Complex{Float64}})
	return get_J(angle.(V),abs.(V),imag.(Y),real.(Y))
end

function get_Jr(θ::Vector{Float64}, V::Vector{Float64}, B::Matrix{Float64}, G::Matrix{Float64})
	n = length(θ)
	Id = diagm(0 => ones(n))
	dI = 1 .- Id

	s = sin.(θ)
	c = cos.(θ)
	vs = V.*s
	vc = V.*c

	Jpt = dI.*(G.*(vs*vc' - vc*vs') - B.*(vc*vc' + vs*vs'))
	Jpt -= diagm(0 => vec(Jpt*ones(n)))
	
	Jqt = -dI.*(G.*(vc*vc' + vs*vs') + B.*(vs*vc' - vc*vs'))
	Jqt -= diagm(0 => vec(Jqt*ones(n)))

	Jpv = dI.*(G.*(vc*c' + vs*s') + B.*(vs*c' - vc*s'))
	Jpv += diagm(0 => (2*V.*diag(G) + Jpv*ones(n)))
	
	Jqv = dI.*(G.*(vs*c' - vc*s') - B.*(vc*c' + vs*s'))
	Jqv += diagm(0 => (-2*V.*diag(B) + Jqv*ones(n)))

	return Jqv - Jqt*pinv(Jpt)*Jpv
end

function get_Jr(V::Vector{Complex{Float64}}, Y::Matrix{Complex{Float64}})
	return get_Jr(angle.(V),abs.(V),imag.(Y),real.(Y))
end


#################################################################





