using Distributions

include("satellite-setup.jl")

G = 1.
dt = .001 # Time resolution of the universe
dt_mpc = .1 # Time resolution of the MPC

r0 = 2.
θ0 = 2π*rand()-π
ω0 = sqrt(G/r0^3)

u0 = [r0*cos(θ0),r0*sin(θ0),-r0*sin(θ0)*ω0,r0*cos(θ0)*ω0]
u = u0 #+ .5*rand(4)

# Target orbit
rt = 3.
ωt = sqrt(G/rt^3)

s = .1
# Parameters error (amplitude,bias)
ξG = (1.,0.)
#ξG = (1.1,0.)
ξr = (1.,0.)
ξθ = (1.,0.)
ξvr = (1.,0.)
ξω = (1.,0.)
ξx = (1.,0.)
#ξx = (1.1,0.)
ξy = (1.,0.)
ξy = ξx
ξvx = (1.,0.)
ξvy = (1.,0.)
ξu = [ξx,ξy,ξvx,ξvy]
ξz = [ξr,ξθ,ξvr,ξω]

figure("Satellite - MPC",(8,8))
PyPlot.plot(0.,0.,"ok")
PyPlot.plot(u[1],u[2],".",color="C0")
PyPlot.plot(rt*cos.(LinRange(-π,π,300)),rt*sin.(LinRange(-π,π,300)),"--k")

# Initial trajectory
T1 = 100.
T = 0.
while T < T1
	global T += 1.
	global u = satellite(u,(0.,0.,G,0.))
	PyPlot.plot(u[1],u[2],".",color="C0")
end

# #=
T2 = 20.

uh = errx(u,ξu)
rh,θh,drh,dθh = radial_coord(uh)
Gh = errx(G,ξG)
ωh = sqrt(Gh/rt^3)
count = 0

while T < T1+T2 && max(abs(rh-rt),abs(dθh-ωt)) > .001
	global T += dt_mpc
	global u

	global count += 1
	if count%10 == 0
		@info "$T/$(T1+T2)"
	end

	global uh = errx(u,ξu)
	global rh,θh,drh,dθh = radial_coord(uh)

	global thr,thθ = thrust_mpc(uh,rt,ωh,Gh,dt_mpc,dt,1.)

	global u = satellite(u,(thr,thθ,G,1.),dt,dt_mpc)
	PyPlot.plot(u[1],u[2],".",color="C1")
end
# =#

T3 = 1000.

while T < T1+T2+T3
	global T += 1.
	global u = satellite(u,(0.,0.,G,1.))
	PyPlot.plot(u[1],u[2],".",color="C2")
end

#title("error factor: $(ξx[1])")
xticks([])
yticks([])



