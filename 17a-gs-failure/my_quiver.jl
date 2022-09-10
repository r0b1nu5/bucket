using PyPlot

function xdot(x::Float64, y::Float64, b::Float64=.91)
	return -3*(b + sin(x)*cos(y))
end

function ydot(x::Float64, y::Float64, b::Float64=.91)
	return -cos(x)*sin(y)
end

function my_quiver(xs::Array{Float64,1}, ys::Array{Float64,1}, scaling::Float64=1., fignum::Int64=111)
	nx = length(xs)
	ny = length(ys)

	U = zeros(nx,ny)
	V = zeros(nx,ny)

	for i in 1:nx
		for j in 1:ny
			x = xs[i]
			y = ys[j]

			U[j,i] = xdot(x,y)
			V[j,i] = ydot(x,y)
		end
	end
	
	figure(fignum)
	PyPlot.quiver(xs,ys,U,V)
end


