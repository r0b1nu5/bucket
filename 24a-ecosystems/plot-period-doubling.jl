using PyPlot, DelimitedFiles

#T = 1_000_000 - 100_000
#δt = 1e-3

#pj3
#T = 19999-1999
#δt = .1

#pj4
#T = 99999 - 1999
#δt = .1

#pj5
T = 19999 - 1999
δt = .1

ωs = (0:T-1)*2π/(T*δt)
F = zeros(0,399)

 #=
nσ = 191
σs = LinRange(2.7,3.7,nσ)

#nσ = 61
#σs = σs[1:61]

for i in 1:nσ
	f = readdlm("data/ex-limit-cycle-04-s$i-fft.csv",',')[[1,],2:200]
	global F = vcat(F,f./maximum(f))
end
# =#

# #=
#nσ = 191
#σs = LinRange(2.7,3.4,nσ)
#ex = "pj2"
##nσ = 9
#σs = [2.,2.32,2.35,2.384,2.417,2.45,2.7,2.776,2.857]
#ex = "pj3"
#ex = "pj4"
#nσ = 20
#σs = [2.0,2.1,2.2,2.3,2.4,2.5,2.6,2.65,2.7,2.75,2.8,2.85,2.9,3.0,3.1,3.2,3.3,3.4,3.5,3.6]

ex = "pj5"
nσ = 5
σs = [2.2,2.4,2.6,2.8,3.0]

for i in 1:nσ
	f = readdlm("data/ex-limit-cycle-"*ex*"-s$i-fft.csv",',')[2:400]'
	global F = vcat(F,f./maximum(f))
end
# =#

figure("Fourier",(12,4))
subplot(2,1,1)
PyPlot.contourf(ωs[2:400],σs,F,50)
for i in 1:nσ
	PyPlot.plot(ωs[[2,400]],[σs[i],σs[i]],"--w",lw=.5)
end
ylabel("σ")
colorbar(label="|fft|")
subplot(2,1,2)
PyPlot.contourf(ωs[2:400],σs,log.(F),50)
for i in 1:nσ
	PyPlot.plot(ωs[[2,400]],[σs[i],σs[i]],"--w",lw=.5)
end
xlabel("ω [rad/s]")
ylabel("σ")
colorbar(label="log|fft|")

figure()
PyPlot.surf(ωs[2:400],σs,F,cmap=get_cmap("viridis"))

figure()
for i in 1:nσ
	n = vec(readdlm("data/ex-limit-cycle-"*ex*"-s$i-N1.csv",','))
	subplot(ceil(Int64,nσ/3),3,i)
	PyPlot.plot((1:length(n))*δt,n,label="σ = $(round(σs[i],digits=2))")
	legend()
	xlabel("t")
	ylabel("N1")
end

#=
PyPlot.matshow(F[nσ:-1:1,2:151])
xticks(0:5:149,round.(ωs[2:5:151],digits=2))
xlabel("ω")
yticks([0,nσ-1],[2.8,2.7])
ylabel("σ")

PyPlot.matshow(log.(F[nσ:-1:1,2:101]))
xticks(0:5:99,round.(ωs[2:5:101],digits=2))
xlabel("ω")
yticks([0,nσ-1],[2.8,2.7])
ylabel("σ")
=#

