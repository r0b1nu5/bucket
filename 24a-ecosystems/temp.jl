N = lv_bunin(N0 + 1e-2*randn(S),A,κ,μ/Si,σ/sqrt(Si),500_000,100_000,1e-3)
surv = vec(1:S)[N[:,end] .> zer0]

figure()
for s in 1:S
	PyPlot.plot(N[s,:])
end
title("np = $np, # survivors = $(length(surv))")





