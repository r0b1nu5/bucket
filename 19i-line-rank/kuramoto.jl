using SparseArrays,Dates

include("L2B.jl")

## INPUT
# L: graph Laplacian
# P: natural frequencies
# th0: initial angles
# store_history: if "true", then returns all angles' time evolution
# verb: displays info along the simulation
function kuramoto(L::Array{Float64,2},P::Array{Float64,1},th0::Array{Float64,1},store_history::Bool=false,verb::Bool=true,max_iter::Int=100000,eps::Float64=1e-6,h::Float64=0.1)
	n = size(L)[1]
	
	B,w = L2B(L)
	W = diagm(0 => w)
		
	th1 = zeros(n)
	th2 = copy(th0)
	
	error = 1000.
	
	iter = 0
	
	ths = Array{Float64,2}(undef,n,0)
	ths = [ths th0]
	dths = Array{Float64,2}(undef,n,0)
	
	while iter < max_iter && error > eps
		iter += 1
		if iter%1000 == 0 || iter == max_iter
			@info "$(now()) -- iter = $iter, err = $error"
		end
		
		th1 = copy(th2)
		
		k1 = P - B*W*sin.(transpose(B)*th1)
		k2 = P - B*W*sin.(transpose(B)*(th1+(h/2)*k1))
		k3 = P - B*W*sin.(transpose(B)*(th1+(h/2)*k2))
		k4 = P - B*W*sin.(transpose(B)*(th1+h*k3))
		
		dth = (k1+2*k2+2*k3+k4)./6
		
		th2 = th1 + h*dth
		
		error = maximum(abs.(dth))
		
		if store_history
			ths = [ths th2]
			dths = [dths dth]
		end
	end
	
	if !store_history
		ths = th2
		dths = dth
	end
	return ths,dths,iter
end


## INPUT
# L: graph Laplacian
# m: vector of inertias
# d: vector of dampings
# P: natural frequencies
# th0: initial angles
# omeg0: initial velocities
# store_history: if "true", then returns all angles' time evolution
# verb: displays info along the simulation
function kuramoto2(L::Array{Float64,2},m::Array{Float64,1},d::Array{Float64,1},P::Array{Float64,1},th0::Array{Float64,1},omeg0::Array{Float64,1},store_history::Bool=false,verb::Bool=true,max_iter::Int=100000,eps::Float64=1e-6,h::Float64=0.025)
	n = size(L)[1]
	
	B,w = L2B(L)
	W = diagm(0 => w)
	M = diagm(0 => m)
	D = diagm(0 => d)
	Mi = inv(M)
	
	th1 = zeros(n)
	th2 = copy(th0)
	omeg1 = zeros(n)
	omeg2 = copy(omeg0)
	
	x1 = [th1;omeg1]
	x2 = [th2;omeg2]
	dx = zeros(2*n)

	error = 1000.
	
	iter = 0
	
	xs = Array{Float64,2}(undef,2*n,0)
	xs = [xs [th0;omeg0]]
	dxs = Array{Float64,2}(undef,2*n,0)
	while iter < max_iter && error > eps
		iter += 1
		if iter%1000 == 0 || iter == max_iter
			@info "$(now()) -- iter = $iter, err = $error"
		end
		
		x1 = copy(x2)
		
		k1 = [(x1[(n+1):(2*n)]);(Mi*P - Mi*B*W*sin.(transpose(B)*x1[1:n]) - Mi*D*x1[(n+1):(2*n)])]
		k2 = [(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]);(Mi*P - Mi*B*W*sin.(transpose(B)*(x1[1:n]+(h/2)*k1[1:n])) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]))]
		k3 = [(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]);(Mi*P - Mi*B*W*sin.(transpose(B)*(x1[1:n]+(h/2)*k2[1:n])) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]))]
		k4 = [(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]);(Mi*P - Mi*B*W*sin.(transpose(B)*(x1[1:n]+h*k3[1:n])) - Mi*D*(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]))]
		
		dx = (k1+2*k2+2*k3+k4)/6
		
		x2 = x1 + h*dx
		
		error = maximum(abs.(dx))
		
		if store_history
			xs = [xs x2]
			dxs = [dxs dx]
		end
	end
	
	if !store_history
		xs = x2
		dxs = dx
	end
	return xs,dxs,iter
end

## INPUT
# L: graph Laplacian
# m: vector of inertias
# d: vector of dampings
# P: natural frequencies
# th0: initial angles
# thref: reference angles to which the deviation are compared !!!!!!!!
# omeg0: initial velocities
# store_history: if "true", then returns all angles' time evolution
# verb: displays info along the simulation
function kuramoto2_P(L::SparseMatrixCSC{Float64,Int},m::Array{Float64,1},d::Array{Float64,1},P::Array{Float64,1},th0::Array{Float64,1},omeg0::Array{Float64,1},thref::Array{Float64,1},verb::Bool=true,max_iter::Int=100000,eps::Float64=1e-6,h::Float64=0.025)
	n = size(L)[1]
	
	B,w = L2B(L)
	W = spdiagm(0 => w)
	M = spdiagm(0 => m)
	D = spdiagm(0 => d)
	Mi = spdiagm(0 => 1 ./ m)
	
	th1 = zeros(n)
	th2 = copy(th0)
	omeg1 = zeros(n)
	omeg2 = copy(omeg0)
	
	x1 = [th1;omeg1]
	x2 = [th2;omeg2]
	dx = zeros(2*n)
	
	error = 1000.
	
	iter = 0
	
	P1 = 0.
	P2 = 0.
	
	while iter < max_iter && error > eps
		iter += 1
		if iter%10000 == 0 || iter == max_iter
			@info "$(now()) -- iter = $iter, err = $error"
		end
		
		x1 = copy(x2)
		
		k1 = [(x1[(n+1):(2*n)]);(Mi*P - Mi*B*W*sin.(transpose(B)*x1[1:n]) - Mi*D*x1[(n+1):(2*n)])]
		k2 = [(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]);(Mi*P - Mi*B*W*sin.(transpose(B)*(x1[1:n]+(h/2)*k1[1:n])) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]))]
		k3 = [(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]);(Mi*P - Mi*B*W*sin.(transpose(B)*(x1[1:n]+(h/2)*k2[1:n])) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]))]
		k4 = [(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]);(Mi*P - Mi*B*W*sin.(transpose(B)*(x1[1:n]+h*k3[1:n])) - Mi*D*(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]))]
		
		dx = (k1+2*k2+2*k3+k4)/6
		
		x2 = x1 + h*dx
		
		error = maximum(abs.(dx))
		
		P1 += sum(h*(x2[1:n]-thref).^2)
		P2 += sum(h*x2[(n+1):(2*n)].^2)
	end
	
	return x2,dx,P1,P2,iter
end

## INPUT
# L: graph Laplacian
# m: vector of inertias
# d: vector of dampings
# P: natural frequencies
# th0: initial angles
# omeg0: initial velocities
# store_history: if "true", then returns all angles' time evolution
# verb: displays info along the simulation
function kuramoto2(L::SparseMatrixCSC{Float64,Int},m::Array{Float64,1},d::Array{Float64,1},P::Array{Float64,1},th0::Array{Float64,1},omeg0::Array{Float64,1},store_history::Bool=false,verb::Bool=true,max_iter::Int=100000,eps::Float64=1e-6,h::Float64=0.025)
	n = size(L)[1]
	
	B,w = L2B(L)
	W = spdiagm(0 => w)
	M = spdiagm(0 => m)
	D = spdiagm(0 => d)
	Mi = spdiagm(0 => 1 ./ m)
	
	th1 = zeros(n)
	th2 = copy(th0)
	omeg1 = zeros(n)
	omeg2 = copy(omeg0)
	
	x1 = [th1;omeg1]
	x2 = [th2;omeg2]
	dx = zeros(2*n)
	
	error = 1000.
	
	iter = 0
	
	xs = Array{Float64,2}(undef,2*n,0)
	xs = [xs [th0;omeg0]]
	dxs = Array{Float64,2}(undef,2*n,0)
	while iter < max_iter && error > eps
		iter += 1
		if iter%1000 == 0 || iter == max_iter
			@info "$(now()) -- iter = $iter, err = $error"
		end
		
		x1 = copy(x2)

#=
@info "$n"
@info "$(size(Mi))"
@info "$(size(P))"
@info "$(size(B))"
@info "$(size(W))"
@info "$(size(D))"
=#
	
		k1 = [(x1[(n+1):(2*n)]);(Mi*P - Mi*B*W*sin.(transpose(B)*x1[1:n]) - Mi*D*x1[(n+1):(2*n)])]
		k2 = [(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]);(Mi*P - Mi*B*W*sin.(transpose(B)*(x1[1:n]+(h/2)*k1[1:n])) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]))]
		k3 = [(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]);(Mi*P - Mi*B*W*sin.(transpose(B)*(x1[1:n]+(h/2)*k2[1:n])) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]))]
		k4 = [(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]);(Mi*P - Mi*B*W*sin.(transpose(B)*(x1[1:n]+h*k3[1:n])) - Mi*D*(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]))]
		
		dx = (k1+2*k2+2*k3+k4)/6
		
		x2 = x1 + h*dx
		
		error = maximum(abs.(dx))
		
		if store_history
			xs = [xs x2]
			dxs = [dxs dx]
		end
	end
	
	if !store_history
		xs = x2
		dxs = dx
	end
	return xs,dxs,iter
end

## INPUT
# B: graph incidence matrix
# w: weights of the lines
# m: vector of inertias
# d: vector of dampings
# P: natural frequencies
# th0: initial angles
# omeg0: initial velocities
# store_history: if "true", then returns all angles' time evolution
# verb: displays info along the simulation
function kuramoto2_incidence(B::SparseMatrixCSC{Float64,Int},w::Array{Float64,1},m::Array{Float64,1},d::Array{Float64,1},P::Array{Float64,1},th0::Array{Float64,1},omeg0::Array{Float64,1},store_history::Bool=false,verb::Bool=true,max_iter::Int=100000,eps::Float64=1e-6,h::Float64=0.025)
	n = size(B)[1]
	
	W = spdiagm(0 => w)
	M = spdiagm(0 => m)
	D = spdiagm(0 => d)
	Mi = spdiagm(0 => 1 ./ m)
	
	th1 = zeros(n)
	th2 = copy(th0)
	omeg1 = zeros(n)
	omeg2 = copy(omeg0)
	
	x1 = [th1;omeg1]
	x2 = [th2;omeg2]
	dx = zeros(2*n)
	
	error = 1000.
	
	iter = 0
	
	xs = Array{Float64,2}(undef,2*n,0)
	xs = [xs [th0;omeg0]]
	dxs = Array{Float64,2}(undef,2*n,0)
	while iter < max_iter && error > eps
		iter += 1
		if iter%1000 == 0 || iter == max_iter
			@info "$(now()) -- iter = $iter, err = $error"
		end
		
		x1 = copy(x2)

#=
@info "$n"
@info "$(size(Mi))"
@info "$(size(P))"
@info "$(size(B))"
@info "$(size(W))"
@info "$(size(D))"
=# 
	
		k1 = [(x1[(n+1):(2*n)]);(Mi*P - Mi*B*W*sin.(transpose(B)*x1[1:n]) - Mi*D*x1[(n+1):(2*n)])]
		k2 = [(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]);(Mi*P - Mi*B*W*sin.(transpose(B)*(x1[1:n]+(h/2)*k1[1:n])) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]))]
		k3 = [(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]);(Mi*P - Mi*B*W*sin.(transpose(B)*(x1[1:n]+(h/2)*k2[1:n])) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]))]
		k4 = [(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]);(Mi*P - Mi*B*W*sin.(transpose(B)*(x1[1:n]+h*k3[1:n])) - Mi*D*(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]))]
		
		dx = (k1+2*k2+2*k3+k4)/6
		
		x2 = x1 + h*dx
		
		error = maximum(abs.(dx))
		
		if store_history
			xs = [xs x2]
			dxs = [dxs dx]
		end
	end
	
	if !store_history
		xs = x2
		dxs = dx
	end
	return xs,dxs,iter
end

function kuramoto2_lin(L::Array{Float64,2},m::Array{Float64,1},d::Array{Float64,1},P::Array{Float64,1},th0::Array{Float64,1},omeg0::Array{Float64,1},store_history::Bool=false,verb::Bool=true,max_iter::Int=100000,eps::Float64=1e-6,h::Float64=0.025)
	n = size(L)[1]
	
	M = diagm(0 => m)
	D = diagm(0 => d)
	Mi = inv(M)
	
	th1 = zeros(n)
	th2 = copy(th0)
	omeg1 = zeros(n)
	omeg2 = copy(omeg0)
	
	x1 = [th1;omeg1]
	x2 = [th2;omeg2]
	dx = zeros(2*n)
	
	error = 1000.
	
	iter = 0
	
	xs = Array{Float64,2}(undef,2*n,0)
	xs = [xs [th0;omeg0]]
	dxs = Array{Float64,2}(undef,2*n,0)
	while iter < max_iter && error > eps
		iter += 1
		if iter%1000 == 0 || iter == max_iter
			@info "$(now()) -- iter = $iter, err = $error"
		end
		
		x1 = copy(x2)
		
		k1 = [(x1[(n+1):(2*n)]);(Mi*P - Mi*L*x1[1:n] - Mi*D*x1[(n+1):(2*n)])]
		k2 = [(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]);(Mi*P - Mi*L*(x1[1:n]+(h/2)*k1[1:n]) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]))]
		k3 = [(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]);(Mi*P - Mi*L*(x1[1:n]+(h/2)*k2[1:n]) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]))]
		k4 = [(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]);(Mi*P - Mi*L*(x1[1:n]+h*k3[1:n]) - Mi*D*(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]))]
		
		dx = (k1+2*k2+2*k3+k4)/6
		
		x2 = x1 + h*dx
		
		error = maximum(abs.(dx))
		
		if store_history
			xs = [xs x2]
			dxs = [dxs dx]
		end
	end
	
	if !store_history
		xs = x2
		dxs = dx
	end
	return xs,dxs
end


function kuramoto2_lin(L::SparseMatrixCSC{Float64,Int},m::Array{Float64,1},d::Array{Float64,1},P::Array{Float64,1},th0::Array{Float64,1},omeg0::Array{Float64,1},store_history::Bool=false,verb::Bool=true,max_iter::Int=100000,eps::Float64=1e-6,h::Float64=0.025)
	n = size(L)[1]
	
	M = spdiagm(0 => m)
	D = spdiagm(0 => d)
	Mi = spdiagm(0 => 1 ./ m)
	
	th1 = zeros(n)
	th2 = copy(th0)
	omeg1 = zeros(n)
	omeg2 = copy(omeg0)
	
	x1 = [th1;omeg1]
	x2 = [th2;omeg2]
	dx = zeros(2*n)
	
	error = 1000.
	
	iter = 0
	
	xs = Array{Float64,2}(undef,2*n,0)
	xs = [xs [th0;omeg0]]
	dxs = Array{Float64,2}(undef,2*n,0)
	while iter < max_iter && error > eps
		iter += 1
		if iter%1000 == 0 || iter == max_iter
			@info "$(now()) -- iter = $iter, err = $error"
		end
		
		x1 = copy(x2)
		
		k1 = [(x1[(n+1):(2*n)]);(Mi*P - Mi*L*x1[1:n] - Mi*D*x1[(n+1):(2*n)])]
		k2 = [(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]);(Mi*P - Mi*L*(x1[1:n]+(h/2)*k1[1:n]) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]))]
		k3 = [(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]);(Mi*P - Mi*L*(x1[1:n]+(h/2)*k2[1:n]) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]))]
		k4 = [(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]);(Mi*P - Mi*L*(x1[1:n]+h*k3[1:n]) - Mi*D*(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]))]
		
		dx = (k1+2*k2+2*k3+k4)/6
		
		x2 = x1 + h*dx
		
		error = maximum(abs.(dx))
		
		if store_history
			xs = [xs x2]
			dxs = [dxs dx]
		end
	end
	
	if !store_history
		xs = x2
		dxs = dx
	end
	return xs,dxs
end

function kuramoto2_lin_P(L::SparseMatrixCSC{Float64,Int},m::Array{Float64,1},d::Array{Float64,1},P::Array{Float64,1},th0::Array{Float64,1},omeg0::Array{Float64,1},thref::Array{Float64,1},verb::Bool=true,max_iter::Int=100000,eps::Float64=1e-6,h::Float64=0.025)
	n = size(L)[1]
	
	M = spdiagm(0 => m)
	D = spdiagm(0 => d)
	Mi = spdiagm(0 => 1 ./ m)
	
	th1 = zeros(n)
	th2 = copy(th0)
	omeg1 = zeros(n)
	omeg2 = copy(omeg0)
	
	x1 = [th1;omeg1]
	x2 = [th2;omeg2]
	dx = zeros(2*n)
	
	error = 1000.
	
	iter = 0
	
	P1 = 0.
	P2 = 0.
	while iter < max_iter && error > eps
		iter += 1
		if iter%1000 == 0 || iter == max_iter
			@info "$(now()) -- iter = $iter, err = $error"
		end
		
		x1 = copy(x2)
		
		k1 = [(x1[(n+1):(2*n)]);(Mi*P - Mi*L*x1[1:n] - Mi*D*x1[(n+1):(2*n)])]
		k2 = [(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]);(Mi*P - Mi*L*(x1[1:n]+(h/2)*k1[1:n]) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]))]
		k3 = [(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]);(Mi*P - Mi*L*(x1[1:n]+(h/2)*k2[1:n]) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]))]
		k4 = [(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]);(Mi*P - Mi*L*(x1[1:n]+h*k3[1:n]) - Mi*D*(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]))]
		
		dx = (k1+2*k2+2*k3+k4)/6
		
		x2 = x1 + h*dx
		
		error = maximum(abs.(dx))
		
		P1 += sum(h*(x2[1:n]-thref).^2)
		P2 += sum(h*x2[(n+1):(2*n)].^2)
	end
	
	return x2,dx,P1,P2,iter
end

# Extracts the max local RoCoF
function rocof_kuramoto2_lin(L::SparseMatrixCSC{Float64,Int},m::Array{Float64,1},d::Array{Float64,1},P::Array{Float64,1},th0::Array{Float64,1},omeg0::Array{Float64,1},verb::Bool=true,max_iter::Int=10,eps::Float64=1e-6,h::Float64=0.025)
	n = size(L)[1]
	
	M = spdiagm(0 => m)
	D = spdiagm(0 => d)
	Mi = spdiagm(0 => 1 ./ m)
	
	A = [spzeros(n,n) spdiagm(0 => ones(n));-Mi*L -Mi*D]
	PP = [zeros(n);Mi*P]
	
	th1 = zeros(n)
	th2 = copy(th0)
	omeg1 = zeros(n)
	omeg2 = copy(omeg0)
	
	x1 = [th1;omeg1]
	x2 = [th2;omeg2]
	dx = zeros(2*n)
	
	error = 1000.
	
	iter = 0
	
	xs = Array{Float64,2}(undef,2*n,0)
	xs = [xs [th0;omeg0]]
	dxs = Array{Float64,2}(undef,2*n,0)
	
	rcf = 0.
	idx = 0
		
	while iter < max_iter && error > eps
		iter += 1
		if verb && (iter%1000 == 0 || iter == max_iter)
			@info "$(now()) -- iter = $iter, err = $error"
		end
		
		x1 = copy(x2)
		
		k1 = PP + A*x1
		k2 = PP + A*(x1 + h/2*k1)
		k3 = PP + A*(x1 + h/2*k2)
		k4 = PP + A*(x1 + h*k3)
		
#=
		k1 = [(x1[(n+1):(2*n)]);(Mi*P - Mi*L*x1[1:n] - Mi*D*x1[(n+1):(2*n)])]
		k2 = [(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]);(Mi*P - Mi*L*(x1[1:n]+(h/2)*k1[1:n]) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]))]
		k3 = [(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]);(Mi*P - Mi*L*(x1[1:n]+(h/2)*k2[1:n]) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]))]
		k4 = [(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]);(Mi*P - Mi*L*(x1[1:n]+h*k3[1:n]) - Mi*D*(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]))]
=#
				
		dx = (k1+2*k2+2*k3+k4)/6
		
		x2 = x1 + h*dx
		
		error = maximum(abs.(dx))
		
		ma,k = findmax(abs.(dx[(n+1):(2*n)])) 
		if ma > abs(rcf)
			rcf = dx[n+k]
			idx = iter
		end
	end
	
	return rcf,idx
end





