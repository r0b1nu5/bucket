using PyPlot

include("plot_ntw.jl")
include("generate_time_series.jl")
include("final.jl")

ntw = "ntw5"
#ntw = "ntw10"
#ntw = "ntw20"
#ntw = "uk"

L = readdlm("data/"*ntw*"_lap_mat.csv",',')
if ntw == "uk"
	L = 10*L
end
m = vec(readdlm("data/"*ntw*"_m.csv",','))
Mi = diagm(0 => 1 ./ m)
d = vec(readdlm("data/"*ntw*"_d.csv",','))

Lm = Mi*L
dm = Mi*d

n = length(d)

a0 = .2
f0 = .009
#f0 = .001
p0 = pi/10

T = 100000
dt = .1

sig = mean(m)*ones(n)

reL = Array{Float64,1}()
rea = Array{Float64,1}()
ref = Array{Float64,1}()
rep = Array{Float64,1}()

for i in 1:n
	@info "$i/$n"

	a = zeros(n)
	a[i] = a0
	f = zeros(n)
	f[i] = f0
	p = zeros(n)
	p[i] = p0

	Xs = generate_forced_time_series(ntw, L, m, d, (a,f,p), T, dt, sig)
#	Xs = load_data(ntw, i)

#	Lh,dh,ah,fh,ph = run_location_large_ntw(Xs, dt, 5)
	Lh,dh,ah,fh,ph = run_location_small_ntw(Xs, dt)
	
	dL = Lh - Lm
	push!(reL, sum(dL.^2)/sum(Lm.^2))
	da = ah[i] - a0
	push!(rea, abs(da)/a0)
	am,id = findmax(ah)
	df = fh[id] - f0
	push!(ref, abs(df)/f0)
	dp = ph[id] - p0
	push!(rep, abs(dp)/p0)

	AA = sortslices([abs.(ah) 1:n],dims=1,rev=true)
	PP = AA[:,1]./sum(AA[:,1])
	id1 = Int(AA[1,2])
	id2 = Int(AA[2,2])
	id3 = Int(AA[3,2])

	figure(666)
	subplot(1,3,3)
	PyPlot.plot([i,i],[0,abs(ah[id1])/sum(abs.(ah))],color="C0")
	PyPlot.plot([i,i],[0,abs(ah[id2])/sum(abs.(ah))],color="C1")
	PyPlot.plot([i,i],[0,abs(ah[id3])/sum(abs.(ah))],color="C2")
	PyPlot.plot(i,abs(ah[i])/sum(abs.(ah)),"o",color="C3",markersize=4.)
end

figure(666)

subplot(1,3,1)
plot_ntw(ntw)

subplot(1,3,2)
PyPlot.semilogy(1:n,reL,"o")
PyPlot.semilogy(n .+ (1:n),rea,"o")
PyPlot.semilogy(2*n .+ (1:n),ref,"o")
PyPlot.semilogy(3*n .+ (1:n),rep,"o")





