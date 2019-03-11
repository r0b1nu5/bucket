include("L2B.jl")

function kuramoto(L::Array{Float64,2},P::Array{Float64,1},th0::Array{Float64,1},verb::Bool=true,max_iter::Int=100000,eps::Float64=1e-6,h::Float64=0.1)
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
			@info("iter = $iter, err = $error")
		end
		
		th1 = copy(th2)
		
		k1 = P - B*W*sin.(transpose(B)*th1)
		k2 = P - B*W*sin.(transpose(B)*(th1+(h/2)*k1))
		k3 = P - B*W*sin.(transpose(B)*(th1+(h/2)*k2))
		k4 = P - B*W*sin.(transpose(B)*(th1+h*k3))
		
		dth = (k1+2*k2+2*k3+k4)./6
		
		th2 = th1 + h*dth
		
		error = maximum(abs.(dth))
		
		ths = [ths th2]
		dths = [dths dth]
	end
	
	return ths,dths
end


function kuramoto2(L::Array{Float64,2},m::Array{Float64,1},d::Array{Float64,1},P::Array{Float64,1},th0::Array{Float64,1},omeg0::Array{Float64,1},verb::Bool=true,max_iter::Int=100000,eps::Float64=1e-6,h::Float64=0.025)
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
	
	error = 1000.
	
	iter = 0
	
	xs = Array{Float64,2}(undef,2*n,0)
	xs = [xs [th0;omeg0]]
	dxs = Array{Float64,2}(undef,2*n,0)
	while iter < max_iter && error > eps
		iter += 1
		if iter%1000 == 0 || iter == max_iter
			@info("iter = $iter, err = $error")
		end
		
		x1 = copy(x2)
		
		k1 = [(x1[(n+1):(2*n)]);(P - Mi*B*W*sin.(transpose(B)*x1[1:n]) - Mi*D*x1[(n+1):(2*n)])]
		k2 = [(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]);(P - Mi*B*W*sin.(transpose(B)*(x1[1:n]+(h/2)*k1[1:n])) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]))]
		k3 = [(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]);(P - Mi*B*W*sin.(transpose(B)*(x1[1:n]+(h/2)*k2[1:n])) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]))]
		k4 = [(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]);(P - Mi*B*W*sin.(transpose(B)*(x1[1:n]+h*k3[1:n])) - Mi*D*(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]))]
		
		dx = (k1+2*k2+2*k3+k4)/6
		
		x2 = x1 + h*dx
		
		error = maximum(abs.(dx))
		
		xs = [xs x2]
		dxs = [dxs dx]
	end
	
	return xs,dxs
end


function kuramoto2_lin(L::Array{Float64,2},m::Array{Float64,1},d::Array{Float64,1},P::Array{Float64,1},th0::Array{Float64,1},omeg0::Array{Float64,1},verb::Bool=true,max_iter::Int=100000,eps::Float64=1e-6,h::Float64=0.025)
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
	
	error = 1000.
	
	iter = 0
	
	xs = Array{Float64,2}(undef,2*n,0)
	xs = [xs [th0;omeg0]]
	dxs = Array{Float64,2}(undef,2*n,0)
	while iter < max_iter && error > eps
		iter += 1
		if iter%1000 == 0 || iter == max_iter
			@info("iter = $iter, err = $error")
		end
		
		x1 = copy(x2)
		
		k1 = [(x1[(n+1):(2*n)]);(P - Mi*L*x1[1:n] - Mi*D*x1[(n+1):(2*n)])]
		k2 = [(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]);(P - Mi*L*(x1[1:n]+(h/2)*k1[1:n]) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k1[(n+1):(2*n)]))]
		k3 = [(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]);(P - Mi*L*(x1[1:n]+(h/2)*k2[1:n]) - Mi*D*(x1[(n+1):(2*n)]+(h/2)*k2[(n+1):(2*n)]))]
		k4 = [(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]);(P - Mi*L*(x1[1:n]+h*k3[1:n]) - Mi*D*(x1[(n+1):(2*n)]+h*k3[(n+1):(2*n)]))]
		
		dx = (k1+2*k2+2*k3+k4)/6
		
		x2 = x1 + h*dx
		
		error = maximum(abs.(dx))
		
		xs = [xs x2]
		dxs = [dxs dx]
	end
	
	return xs,dxs
end



