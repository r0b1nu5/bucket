using Statistics

function kuramoto_R2(w::Array{Float64,1}, x0::Array{Float64,1}, y0::Array{Float64,1}, K::Float64=1., h::Float64=.1, max_iter::Int64=1000, eps::Float64=1e-6)
	n = length(x0)
	
	x1 = copy(x0)
	y1 = copy(y0)
	x2 = copy(x0)
	y2 = copy(y0)
	
	err = 1000.
	it = 0
	
	xs = Array{Float64,2}(undef,n,0)
	xs = [xs x0]
	ys = Array{Float64,2}(undef,n,0)
	ys = [ys y0]
	
	while err > eps && it < max_iter
		it += 1
		
		x1 = copy(x2)
		y1 = copy(y2)

		kx1 = fx(x1,y1,w,K)
		ky1 = fy(x1,y1,w,K)
		kx2 = fx(x1+h/2*kx1,y1+h/2*ky1,w,K)
		ky2 = fy(x1+h/2*kx1,y1+h/2*ky1,w,K)
		kx3 = fx(x1+h/2*kx2,y1+h/2*ky2,w,K)
		ky3 = fy(x1+h/2*kx2,y1+h/2*ky2,w,K)
		kx4 = fx(x1+h*kx3,y1+h*ky3,w,K)
		ky4 = fy(x1+h*kx3,y1+h*ky3,w,K)
		
		dx = (kx1 + 2*kx2 + 2*kx3 + kx4)/6
		dy = (ky1 + 2*ky2 + 2*ky3 + ky4)/6
		
		x2 = x1 + h*dx
		y2 = y1 + h*dy
		
		xs = [xs x2]
		ys = [ys y2]
		
		err = max(maximum(abs.(dx)),maximum(abs.(dy)))
	end
	
	return xs,ys
end


function fx(x::Array{Float64,1}, y::Array{Float64,1}, w::Array{Float64,1}, K::Float64=1.)
	r = sqrt.(x.^2 + y.^2)
	x0 = mean(x./r)
	y0 = mean(y./r)
	
	dx = K*x0*y.^2 ./r - K*y0*x.*y./r - w.*y
	
	return dx
end
		

function fy(x::Array{Float64,1}, y::Array{Float64,1}, w::Array{Float64,1}, K::Float64=1.)
	r = sqrt.(x.^2 + y.^2)
	x0 = mean(x./r)
	y0 = mean(y./r)
	
	dy = K*y0*x.^2 ./r - K*x0*x.*y./r + w.*x
	
	return dy
end
