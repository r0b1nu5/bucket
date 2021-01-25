using PyPlot, Statistics

include("ksakaguchi.jl")

phi1 = -.6
phi2 = 1.05
phis = LinRange(phi1,phi2,300)

a = .5
g1 = -1.
g2 = 1.

om = 2*rand(3)
om .-= mean(om)

V = [1,2,3]
E = [(1,3),(2,3),(3,1),(3,2)]

L = [1. 0. -1.;0. 1. -1.;-1. -1. 2.]

n = 3
m = 2

ps = Array{Float64,1}()
fs = Array{Float64,2}(undef,2*m,0)

best_phi = 0.
best_dis = 1000.
best_om2 = 0.

for phi in phis
	global f = zeros(2*m)	
	global om2 = copy(om)

	for i in 1:n-1
		if om2[i]-g2 <= phi <= om2[i]-g1
			f[i] = om2[i] - phi
			f[i+m] = sin(-asin(f[i])-2*a)

			global e = E[i]

			for j in 1:n
				if e[2] == j
					global om2[j] -= f[i+m]
				end
			end
		else
			f = zeros(2*m)
			break
		end
	end
	
	push!(ps,phi)
	global fs = [fs f]

	dis = abs(om2[n] - phi)
	if dis < best_dis
		global best_phi = copy(phi)
		global best_om2 = copy(om2[n])
		global best_dis = copy(dis)
	end
end

ths,dhs,err,iter = ksakaguchi(L,om,zeros(n),a)
ff = Array{Float64,1}()
for e in E
	push!(ff,sin(ths[e[1],end]-ths[e[2],end]-a))
end

PyPlot.plot(ps,ps,"--")
for i in 1:2*m
	PyPlot.plot(ps,fs[i,:],color="C$(i)")
end
#PyPlot.plot(ps,fs[end,:])

if iter < 100000
	d = mean(dhs[:,end])
	for i in 1:2*m
		PyPlot.plot(d,ff[i],"o",color="C$(i)")
	end
#	PyPlot.plot(d,ff[end],"o")
end

PyPlot.plot([best_phi,best_phi],[-1.,1.],"--k")
PyPlot.plot([phi1,phi2],[best_om2,best_om2],"--k")






