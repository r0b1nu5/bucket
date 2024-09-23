using PyPlot

nσ = 20
σs = LinRange(2.7,2.8,nσ)

T = 2_000_000 - 100_000
δt = 1e-3
ωs = (0:T-1)*2π/(T*δt)
F = zeros(0,199)

for i in 1:nσ
	global F = vcat(F,readdlm("data/ex-limit-cycle-04-s$i-fft.csv",',')[[1,],2:200])
end

figure()
subplot(2,1,1)
PyPlot.contourf(ωs[2:200],σs,F)
subplot(2,1,2)
PyPlot.contourf(ωs[2:200],σs,log.(F))

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

