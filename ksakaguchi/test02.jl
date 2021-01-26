using PyPlot, Statistics

include("ksakaguchi.jl")

phi1 = .3
phi2 = .85
phis = LinRange(phi1,phi2,300)

a = .4
g1 = -1.
g2 = sin(pi/2-2*a)

 #=
V = [1,2,3]
E = [(1,3),(2,3),(3,1),(3,2)]
L = [1. 0. -1.;0. 1. -1.;-1. -1. 2.]
n = 3
m = 2
# =#
# #=
V = [1,2,3,4]
E = [(1,3),(2,3),(3,4),(3,1),(3,2),(4,3)]
L = [1. 0. -1. 0.;0. 1. -1. 0.;-1. -1. 3. -1.;0. 0. -1. 1.]
n = 4
m = 3
# =#
 #=
om0 = 2.
om = om0*rand(n)
om .-= mean(om)
# =#

ps = Array{Float64,1}()
fs = Array{Float64,2}(undef,2*m,0)
om2s = Array{Float64,2}(undef,n,0)

best_phi = 0.
best_dis = 1000.
best_om2 = 0.

for phi in phis
	global f = zeros(2*m)	
	global om2 = copy(om)
	local ok = true

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
			om2 = copy(om)
			ok = false
			break
		end
	end
	
	push!(ps,phi)
	global fs = [fs f]
	global om2s = [om2s om2]


	local dis = abs(om2[n] - phi)
	if dis < best_dis && ok
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

figure()

subplot(1,2,1)
for i in 1:m
	PyPlot.plot(ps,fs[i,:],color="C$(i-1)")
	PyPlot.plot(ps,fs[i+m,:],"--",color="C$(i-1)")

	if iter < 100000
		global d = mean(dhs[:,end])
		PyPlot.plot(d,ff[i],"o",color="C$(i-1)")
		PyPlot.plot(d,ff[i+m],"s",color="C$(i-1)")
	end
end
xlabel("φ (flow parameter)")
ylabel("f_ij")

lb = min(minimum(ps),minimum(om2s[n,:]))
ub = max(maximum(ps),maximum(om2s[n,:]))

subplot(1,2,2)
PyPlot.plot(ps,ps)
PyPlot.plot(ps,om2s[n,:],".")
PyPlot.plot([best_phi,best_phi],[lb,ub],"--k")
PyPlot.plot([phi1,phi2],[best_om2,best_om2],"--k")
xlabel("φ (flow parameter)")
ylabel("ω_n^{(n)}")





