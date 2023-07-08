using LinearAlgebra, PyPlot

# u=(x,y,dx,dy): cartesian position and velocity
# p=(r,θ,G,maxth):
# 	r: amplitude of the thrust
# 	θ: angle of the thrust
# 	G: effective gravitational constant
# 	maxth: maximal thrust
# t:time
function Fg(u,p,t=0.)
	fx = p[4]*p[1]*cos(p[2])
	fy = p[4]*p[1]*sin(p[2])
	
	r = sqrt(u[1]^2+u[2]^2)
	n = -[u[1],u[2]]./r

	du = zeros(4)
	du[1] = u[3]
	du[2] = u[4]
	du[3] = p[3]*n[1]/r^2 + fx
	du[4] = p[3]*n[2]/r^2 + fy

	return du
end

# u0=(x0,y0,dx0,dy0): initial position and velocity
# p=(r,θ,G,maxth): parameters
# dt: time step
# T: simulation time
function satellite(u0,p,dt=.01,T=1.)
	N = T/dt
	u = u0

	for t in 1:N
		k1 = Fg(u,p,t*dt)
		k2 = Fg(u+dt*k1/2,p,(t+.5)*dt)
		k3 = Fg(u+dt*k2/2,p,(t+.5)*dt)
		k4 = Fg(u+dt*k3,p,(t+1)*dt)

		du = (k1+2*k2+2*k3+k4)/6

		u += dt*du
	end

	return u
end

function err_u(u,ξu)
	uh = [ξu[1][1]*u[1]+ξu[1][2],
	      ξu[2][1]*u[2]+ξu[2][2],
	      ξu[3][1]*u[3]+ξu[3][2],
	      ξu[4][1]*u[4]+ξu[4][2]]

	return uh
end

# Tries to determine the optimal thrust to reach ut from uh.
function thrust_mpc(uh,rt,ωt,Gh,dt_mpc,dt,thrmax=1.)
	rh,θh,drh,dθh = radial_coord(uh)
	
	thr = [0.,θh] # amplitude and angle of the thrust initialized at 0 and angle of the position vector. 
	u = satellite(uh,(thr[1],thr[2],Gh,1.),dt,10*dt_mpc)
	r,θ,dr,dθ = radial_coord(u)
	
	δr0 = rt - rh
	δr1 = rt - r
	δω0 = ωt - dθh
	δω1 = ωt - dθ

	counter = 0
	incr = .01*thrmax
	tr = Inf
	tθ = Inf
	trnew = 0.
	tθnew = θh

	while (δr0*δr1 > 0 || δω0*δω1 > 0) && (abs(tr-trnew) > incr/2 || abs(tθ-tθnew) > incr/2) && (counter < 1e5)
		counter += 1
#		if counter%100 == 0
#			@info "counter = $counter, tr = $tr, tθ = $tθ"
#			@info "$(δr0*δr1 > 0), $(δω0*δω1 > 0), $(abs(tr-trnew) > incr/2), $(abs(tθ-tθnew) > incr/2)"
#			@info "$((δr0*δr1 > 0 || δω0*δω1 > 0) && (abs(tr-trnew) > incr/2 || abs(tθ-tθnew) > incr/2) && counter < 1e5)"
#		end

		tr = trnew
		tθ = tθnew
		
		if δr0*δr1 > 0
			trnew = clamp(tr + incr*sign(δr0),-thrmax,thrmax)
		end
		if δω0*δω1 > 0
			tθnew = clamp(tθ + incr*sign(δω0)*sign(δr0),θh-π/2,θh+π/2)
		end

		u = satellite(uh,[trnew;tθnew;Gh;1.],dt,10*dt_mpc)
		r,θ,dr,dθ = radial_coord(u)
		
		δr0 = rt - rh
		δr1 = rt - r
		δω0 = ωt - dθh
		δω1 = ωt - dθ
	end

	tr = trnew
	tθ = tθnew

	if counter >= 1e5
		@warn "Max. counter reached"
	end

	return tr,tθ
end


function radial_coord(u)
    x, y, vx, vy  = u[1], u[2], u[3], u[4]
    r = sqrt(x^2 + y^2)
    θ = atan(y, x)
    drdt = (x * vx + y * vy) / sqrt(x^2 + y^2)
    dθdt = (vy * x - y * vx) / (y^2 + x^2)
    return [r; θ; drdt; dθdt]
end

function cart_coord(v)
    r, θ, drdt, dθdt  = v[1], v[2], v[3], v[4]
    x = r * cos(θ)
    y = r * sin(θ)
    vx = drdt * cos(θ) - r * sin(θ) * dθdt
    vy = drdt * sin(θ) + r * cos(θ) * dθdt
    return [x; y; vx; vy]
end





