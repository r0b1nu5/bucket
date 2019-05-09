using LinearAlgebra

include("L2B.jl")

# B: incidence matrix of the interaction graph!!!
# w: edge weights
function NR_kuramoto(B::Array{Float64,2},w::Array{Float64,1},omega::Array{Float64,1},eps::Float64=1e-8,max_iter::Int64=50)
	W = diagm(0 => w)
	Bt = transpose(B)
	L = B*W*Bt
	A = L.*(diagm(0 => ones(length(omega))) .- 1)
	
	th = zeros(length(omega))
	th .-= th[end]
	f = omega - B*W*sin.(Bt*th)
	
	iter = 0
	err = maximum(abs.(f))
	
	while err > eps && iter < max_iter
		iter += 1
		
		C = cos.(th)
		S = sin.(th)
		J = A.*(C*transpose(C)+S*transpose(S))
		J -= diagm(0 => vec(sum(J,dims=2)))
		
# We work with the system reduced by one vertex to remove the invariance by a constant phase shift. Otherwise, the Jacobian is not invertible.
		thr = th[1:end-1]
		Jr = J[1:end-1,1:end-1]
		fr = f[1:end-1]
		
		thr = thr - inv(Jr)*fr
		th = [thr;0]
		
		f = omega - B*W*sin.(Bt*th)
		err = maximum(abs.(f))
		
		@info "iter = $iter, error = $err"
	end
	
	if iter == max_iter
		@info "No convergence..."
	end
	
	return th,iter
end

function NR_kuramoto(B::SparseMatrixCSC{Float64,Int64},w::Array{Float64,1},omega::Array{Float64,1},eps::Float64=1e-8,max_iter::Int64=50)
	W = spdiagm(0 => w)
	Bt = transpose(B)
	L = B*W*Bt
	A = L.*(spdiagm(0 => ones(length(omega))) .- 1)
	
	th = zeros(length(omega))
	th .-= th[end]
	f = omega - B*W*sin.(Bt*th)
	
	iter = 0
	err = maximum(abs.(f))
	
	while err > eps && iter < max_iter
		iter += 1
		
		C = cos.(th)
		S = sin.(th)
		J = A.*(C*transpose(C)+S*transpose(S))
		J -= diagm(0 => vec(sum(J,dims=2)))
		
# We work with the system reduced by one vertex to remove the invariance by a constant phase shift. Otherwise, the Jacobian is not invertible.
		thr = th[1:end-1]
		Jr = J[1:end-1,1:end-1]
		fr = f[1:end-1]
		
		thr = thr - inv(Jr)*fr
		th = [thr;0]
		
		
		f = omega - B*W*sin.(Bt*th)
		err = maximum(abs.(f))
		
		@info "iter = $iter, error = $err"
	end
	
	if iter == max_iter
		@info "No convergence..."
	end
	
	return th,iter
end

function NR_kuramoto(B::Array{Float64,2},w::Array{Float64,1},omega::Array{Float64,1},th0::Array{Float64,1},eps::Float64=1e-8,max_iter::Int64=50)
	W = diagm(0 => w)
	Bt = transpose(B)
	L = B*W*Bt
	A = L.*(diagm(0 => ones(length(omega))) .- 1)
	
	th = th0 .- th0[end]
	f = omega - B*W*sin.(Bt*th)
	
	iter = 0
	err = maximum(abs.(f))
	
	while err > eps && iter < max_iter
		iter += 1
		
		C = cos.(th)
		S = sin.(th)
		J = A.*(C*transpose(C)+S*transpose(S))
		J -= diagm(0 => vec(sum(J,dims=2)))
		
# We work with the system reduced by one vertex to remove the invariance by a constant phase shift. Otherwise, the Jacobian is not invertible.
		thr = th[1:end-1]
		Jr = J[1:end-1,1:end-1]
		fr = f[1:end-1]
		
		thr = thr - inv(Jr)*fr
		th = [thr;0]
		
		
		f = omega - B*W*sin.(Bt*th)
		err = maximum(abs.(f))
		
		@info "iter = $iter, error = $err"
	end
	
	if iter == max_iter
		@info "No convergence..."
	end
	
	return th,iter
end

function NR_kuramoto(B::SparseMatrixCSC{Float64,Int64},w::Array{Float64,1},omega::Array{Float64,1},th0::Array{Float64,1},eps::Float64=1e-8,max_iter::Int64=50)
	W = spdiagm(0 => w)
	Bt = transpose(B)
	L = B*W*Bt
	A = L.*(spdiagm(0 => ones(length(th0))) .- 1)
	
	th = th0 .- th0[end]
	f = omega - B*W*sin.(Bt*th)
	
	iter = 0
	err = maximum(abs.(f))
	
	while err > eps && iter < max_iter
		iter += 1
		
		C = cos.(th)
		S = sin.(th)
		J = A.*(C*transpose(C)+S*transpose(S))
		J -= diagm(0 => vec(sum(J,dims=2)))
		
# We work with the system reduced by one vertex to remove the invariance by a constant phase shift. Otherwise, the Jacobian is not invertible.
		thr = th[1:end-1]
		Jr = J[1:end-1,1:end-1]
		fr = f[1:end-1]
		
		thr = thr - inv(Jr)*fr
		th = [thr;0]
		
		
		f = omega - B*W*sin.(Bt*th)
		err = maximum(abs.(f))
		
		@info "iter = $iter, error = $err"
	end
	
	if iter == max_iter
		@info "No convergence..."
	end
	
	return th,iter
end

function isstable(L::Array{Float64,2},th::Array{Float64,1})
	A = L.*(diagm(0 => ones(length(th))) .- 1)
	C = cos.(th)
	S = sin.(th)
	J = A.*(C*transpose(C)+S*transpose(S))
	J -= diagm(0 => vec(sum(J,dims=2)))
	
	ls = eigvals(J)
	
	t = false
	if maximum(ls) < 1e-8
		t = true
	end
	
	return t
end

function isstable(L::SparseMatrixCSC{Float64,Int64},th::Array{Float64,1})
	A = L.*(spdiagm(0 => ones(length(th))) .- 1)
	C = cos.(th)
	S = sin.(th)
	J = A.*(C*transpose(C)+S*transpose(S))
	J -= spdiagm(0 => vec(sum(J,dims=2)))
	
	ls = eigvals(Array(J))
	
	t = false
	if maximum(ls) < 1e-8
		t = true
	end
	
	return t
end
	







