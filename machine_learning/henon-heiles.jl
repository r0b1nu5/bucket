using PyPlot, DelimitedFiles

# State of the system in encoded as xy = (x,\dot{x},y,\dot{y})

function hh(xy0::Array{Float64,1},do_the_plot::Bool=false,n_iter::Int64=1000,h::Float64=.01,lambda::Float64=1.)
	x0,xd0,y0,yd0 = xy0

	E = E_hh(xy0,lambda)
	@info "Energy level: E = $E"

	iter = 0
	c = 0
	xy1 = copy(xy0)
	xy2 = copy(xy0)
	xys = Array{Float64,2}(undef,4,0)
	xys = [xys xy0]
	
	while iter < n_iter && abs(E_hh(xy1,lambda)-E) < 1e-8
		iter += 1
		if iter%1000 == 0
			@info "$iter/$n_iter"
		end
		
		xy1 = copy(xy2)

		k1 = f(xy1,lambda)
		k2 = f(xy1 + h/2*k1,lambda)
		k3 = f(xy1 + h/2*k2,lambda)
		k4 = f(xy1 + h*k3,lambda)

		dxy = (k1+2*k2+2*k3+k4)/6
		xy2 = xy1 + h*dxy

		xys = [xys xy2]

		if iter%10000 == 0
			c += 1
			writedlm("data1/xys$c.csv",xys,',')
			xys = Array{Float64,2}(undef,4,0)
			xys = [xys xy2]
		end
	end

	XYs = Array{Float64,2}(undef,4,0)
	for i in 1:c
		xy = readdlm("data1/xys$i.csv",',')
		XYs = [XYs xy[:,1:end-1]]
	end
	XYs = [XYs xys]
	xys = XYs

	if abs(E_hh(xys[:,end],lambda) - E) > 1e-8
		@info "Diverged (energy not conserved)."
		if do_the_plot
			PyPlot.plot(xys[1,1],xys[3,1],"o",color="C4")
		end
	else
		if do_the_plot
			PyPlot.plot(xys[1,:],xys[3,:],color="C0")
			PyPlot.plot(xys[1,1],xys[3,1],"o",color="C3")
			PyPlot.plot(xys[1,end],xys[3,end],"o",color="C2")
		end
	end

	return xys, E
end


function V_hh(xy::Array{Float64,1},lambda::Float64=1.)
	x,xd,y,yd = xy

	return .5*(x^2 + y^2) + lambda*(x^2*y - y^3/3)
end

function E_hh(xy::Array{Float64,1},lambda::Float64=1.)
	x,xd,y,yd = xy

	return V_hh(xy,lambda) + .5*(xd^2 + yd^2)
end

function f(xy::Array{Float64,1},lambda::Float64=1.)
	x,xd,y,yd = xy

	return [xd,
		-x - 2*lambda*x*y,
		yd,
		-y - lambda*(x^2 - y^2)]
end


