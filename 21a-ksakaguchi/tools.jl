using Statistics, LinearAlgebra, Dates, SparseArrays

function winding(θ::Array{Float64,1}, Σ::Array{Int64,1})
	if length(θ) < 1
		q = []
	else
		θ1 = θ[Σ]
		θ2 = θ[[Σ[2:end];Σ[1]]]

		dθ = θ1 - θ2

		q = round(Int,sum(mod.(dθ .+ π,2π) .- π)/(2π))
	end

	return q
end

# Computes the winding vector for cycles defined by the cycles in C.
function winding(θ::Array{Float64,1}, Σ::Array{Array{Int64,1},1})
	if length(Σ) < 1 || length(θ) < 1
		qs = []
	else
		qs = [winding(θ,Σ[i]) for i in 1:length(Σ)]
	end

	return qs
end

# Computes the max (in abs value) angle difference, using the incidence matrix.
function cohesiveness_inc(θ::Array{Float64,1}, B::Array{Float64,2})
	return maximum(abs.(mod.(transpose(B)*θ .+ π,2π) .- π))
end

# Computes the max (in abs value) angle difference, using the adjacency/Laplacian matrix.
function cohesiveness_adj(θ::Array{Float64,1}, A::Array{Float64,2})
	n = length(θ)
	
	dθ = (θ*ones(1,n) - ones(n)*θ').*(abs.(A .* (1 .- diagm(0 => ones(n)))) .> 1e-8)

	return maximum(dθ)
end

# Draws a random element of the winding cell of u
# u: Winding vector
# B: Incidence matrix (n x m)
# T: Indices of the edges composing a spanning tree
# C: cycle-edge incidence matrix (see Jafarpour et al. (2020))
# TODO Are we really sampling the whole winding cell???
function sample_winding_cell(u::Array{Int64,1}, B::Array{Float64,2}, T::Array{Int64,1}, C::Array{Float64,2}, max_attempt::Int64=1000)
	n,m = size(B)

	Bt = transpose(B)
	Ttd = pinv(transpose(B[:,T]))
	Cd = pinv(C)

	c = 0
	test = true
	θ = Array{Float64,1}()

	while c < max_attempt && test
		c += 1

		x = 2π*rand(n) .- π
		x .-= mean(x)

		y = Bt*x + 2π*Cd*u

		if norm(y,Inf) < π
			θ = Ttd*y[T]
			θ = θ .- mean(θ)
			test = false
		end
	end

	@info "Number of attempts: $c"
	if c == max_attempt
		@info "WARNING: No sample found."
	end

	return θ
end

function jacobian(L::Array{Float64,2}, θ::Array{Float64,1}, α::Float64)
	n = length(θ)

	A = -L.*(1 .- diagm(0 => ones(n)))

	dθ = θ*ones(1,n) - ones(n)*θ'

	J = A.*cos.(dθ .- α)
	J = J - diagm(0 => vec(sum(J,dims=2)))

	return J
end


function cyqle(n::Int64, q::Int64=1)
	A = zeros(n,n)
	for i in 1:q
		A += diagm(i => ones(n-i)) + diagm(-i => ones(n-i)) + diagm(n-i => ones(i)) + diagm(i-n => ones(i))
	end

	D = diagm(0 => vec(sum(A,dims=2)))

	L = D - A

	return L
end

function gershgorin(M::Union{Array{Float64,2},Array{Complex{Float64},2}})
	c = diag(M)

	n = length(c)

	M0 = (1 .- diagm(0 => ones(n))).*M

	rh = vec(sum(abs.(M0),dims=2))
	rv = vec(sum(abs.(M0),dims=1))

	ls = eigvals(M)

	t = LinRange(0,2pi,100)
	x1 = minimum(real.(c)-rh)-1
	x2 = maximum(real.(c)+rh)+1
	y1 = minimum(imag.(c)-rh)-1
	y2 = maximum(imag.(c)+rh)+1
	x3 = minimum(real.(c)-rv)-1
	x4 = maximum(real.(c)+rv)+1
	y3 = minimum(imag.(c)-rv)-1
	y4 = maximum(imag.(c)+rv)+1

	figure()
	subplot(1,2,1)
	PyPlot.plot([x1,x2],[0.,0.],"--k")
	PyPlot.plot([0.,0.],[y1,y2],"--k")
	axis([x1,x2,y1,y2])
	subplot(1,2,2)
	PyPlot.plot([x3,x4],[0.,0.],"--k")
	PyPlot.plot([0.,0.],[y3,y4],"--k")
	axis([x3,x4,y3,y4])

	for i in 1:n
		subplot(1,2,1)
		PyPlot.plot(rh[i]*cos.(t) .+ real(c[i]),rh[i]*sin.(t) .+ imag(c[i]),label="i = $i")
		PyPlot.plot(real(ls[i]),imag(ls[i]),"ok")
		subplot(1,2,2)
		PyPlot.plot(rv[i]*cos.(t) .+ real(c[i]),rv[i]*sin.(t) .+ imag(c[i]))
		PyPlot.plot(real(ls[i]),imag(ls[i]),"ok")
	end

	subplot(1,2,1)
	title("Row Gershgorin circles")
	legend()
	subplot(1,2,2)
	title("Column Gershgoring circles")

	return c,rh,rv
end

function dcc(x::Float64)
	return mod(x + π,2π) - π
end

function dcc(x::Array{Float64,1})
	return [dcc(x[i]) for i in 1:length(x)]
end

function dcc(x::Array{Float64,2})
	d = Array{Float64,2}(undef,size(x)[1],0)
	for j in 1:size(x)[2]
		d = [d dcc(x[:,j])]
	end
	
	return d
end

function dcc2(x::Float64)
	return mod(x + π/2,2π) - π/2
end

function dcc2(x::Array{Float64,1})
	return [dcc2(x[i]) for i in 1:length(x)]
end

function dcc2(x::Array{Float64,2})
	d = Array{Float64,2}(undef,size(x)[1],0)
	for j in 1:size(x)[2]
		d = [d dcc2(x[:,j])]
	end
	
	return d
end

function winding2(θ::Array{Float64,1}, Σ::Array{Int64,1})
	if length(θ) < 1
		q = []
	else
		θ1 = θ[Σ]
		θ2 = θ[[Σ[2:end];Σ[1]]]

		dθ = θ1 - θ2

		q = round(Int,sum(dcc2(dθ))/(2π))
	end

	return q
end

function L2B(L::Array{Float64,2})
	n = size(L)[1]
	
	B = Array{Float64,2}(undef,n,0)
	w = Array{Float64,1}(undef,0)
	for i in 1:n-1
		for j in i+1:n
			if L[i,j] != 0.0
				ed = zeros(n)
				ed[i] = 1.0
				ed[j] = -1.0
				B = [B ed]
				push!(w,-L[i,j])
			end
		end
	end
	
	return B,w
end

function L2B(L::SparseMatrixCSC{Float64,Int})
	n = size(L)[1]
	
	B = Array{Float64,2}(undef,n,0)
	w = Array{Float64,1}(undef,0)
	Bt = Array{Float64,2}(undef,0,n)
	for i in 1:n-1
		for j in i+1:n
			if L[i,j] != 0.0
				ed = zeros(n)
				edt = zeros(1,n)
				ed[i] = 1.0
				edt[1,i] = 1.0
				ed[j] = -1.0
				edt[1,j] = -1.0
				B = [B ed]
				Bt = [Bt;edt]
				push!(w,-L[i,j])
			end
		end
	end
	
	return sparse(B),w,sparse(Bt)
end

function L2B_ij(L::Array{Float64,2}, i0::Int64, j0::Int64)
	n = size(L)[1]

	if i0 > j0
		k = i0
		i0 = copy(j0)
		j0 = copy(i0)
	end

	B = Array{Float64,2}(undef,n,0)
	w = Array{Float64,1}(undef,0)
	c = 0
	eij = 0
	for i in 1:n-1
		for j in i+1:n
			if L[i,j] != 0.0 && i == i0 && j == j0
				c += 1
				eij = copy(c)
				ed = zeros(n)
				ed[i] = 1.0
				ed[j] = -1.0
				B = [B ed]
				push!(w,-L[i,j])
			elseif L[i,j] != 0.0
				c += 1
				ed = zeros(n)
				ed[i] = 1.0
				ed[j] = -1.0
				B = [B ed]
				push!(w,-L[i,j])
			end
		end
	end

	return B,w,eij
end

function L2B_bidir(L::Array{Float64,2})
	b,w = L2B(L)

	B1 = b.*(b .> 0)
	B2 = -b.*(b .< 0)

	Bout = [B1 B2]
	Bin = [B2 B1]

	B = Bout - Bin

	return B,Bout,Bin
end


# Computes the out-incidence matrix of the bidirected counterpart of an undirected graph, based on the incidence matrix of the latter.

function B2Bout(B::Array{Float64,2})
	B1 = B.*(B .> 0.)
	B2 = -B.*(B .< 0.)

	return [B1 B2]
end

# Same for the in-incidence matrix.

function B2Bin(B::Array{Float64,2})
	B1 = B.*(B .> 0.)
	B2 = -B.*(B .< 0.)

	return [B1 B2]
end

function proj_mtx(B::Array{Float64,2}, Bout::Array{Float64,2}, D::Array{Float64,2})
	n,m = size(B)

	Id = diagm(0 => ones(m))

	PD = Id - D*B'*pinv(Bout*D*B')*Bout

	Q = 2*PD'*PD + Id - PD - PD'

	return PD,Q
end


function rand_Q(n::Int64)
	R = 2*rand(n,n) .- 1
	while minimum(abs.(eigvals(R))) < 1e-8
		R = 2*rand(n,n) .- 1
	end

	return R'*R
end


function find_ortho_vec(A::Array{Float64,2})
	n,m = size(A)

	ma = -1.
	v = zeros(n)

	while ma < 1e-8
		v = 2*rand(n) .- 1

		for i in 1:m
			u = A[:,i]
			u /= norm(u)
			v -= dot(v,u)*u
		end

		v /= norm(v)
		ma = maximum(abs.(v))
	end

	return v
end

function gen_R(n::Int64)
	R = zeros(n-1,n)
	for i in 1:n-1
		R[i,1:i] = ones(1,i)./sqrt(i+i^2)
		R[i,i+1] = -i/sqrt(i+i^2)
	end

	return R
end

