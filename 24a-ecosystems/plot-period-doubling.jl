using PyPlot, DelimitedFiles

nσ = 191
σs = LinRange(2.7,3.7,nσ)

T = 2_000_000 - 100_000
δt = 1e-3
ωs = (0:T-1)*2π/(T*δt)
F = zeros(0,199)

#nσ = 61
#σs = σs[1:61]

for i in 1:nσ
	f = readdlm("data/ex-limit-cycle-04-s$i-fft.csv",',')[[1,],2:200]
	global F = vcat(F,f./maximum(f))
end

figure("Fourier",(12,4))
subplot(2,1,1)
PyPlot.contourf(ωs[2:200],σs,F,50)
ylabel("σ")
colorbar(label="|fft|")
subplot(2,1,2)
PyPlot.contourf(ωs[2:200],σs,log.(F),50)
xlabel("ω [rad/s]")
ylabel("σ")
colorbar(label="log|fft|")

figure()
PyPlot.surf(ωs[2:200],σs,F,cmap=get_cmap("viridis"))
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

