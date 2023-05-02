using PyPlot, DelimitedFiles

n = 20
ks = 1:40
τ = .1
T1 = 100
T2 = 230

function rescale(x::Union{Matrix{Float64}, Vector{Float64}})
	return (abs.(x) .- minimum(abs.(x)))./(maximum(abs.(x)) - minimum(abs.(x)))
end

L100l1 = vec(readdlm("test-nonagnostic/l1/Agnostic_l1_100_freq",','))
γ100l1 = vec(readdlm("test-nonagnostic/l1/Agnostic_l1_100_gamma",','))
L230l1 = vec(readdlm("test-nonagnostic/l1/Agnostic_l1_230_freq",','))
γ230l1 = vec(readdlm("test-nonagnostic/l1/Agnostic_l1_230_gamma",','))

L100l1_Aa = vec(readdlm("test-nonagnostic/l1/l1_100_freq",','))
γ100l1_Aa = vec(readdlm("test-nonagnostic/l1/l1_100_gamma",','))
L230l1_Aa = vec(readdlm("test-nonagnostic/l1/l1_230_freq",','))
γ230l1_Aa = vec(readdlm("test-nonagnostic/l1/l1_230_gamma",','))

figure("SALO-relax")

subplot(2,4,1)
PyPlot.plot(ks./(T1*τ),rescale(L100l1),color="C0")
xlabel("f")
title("rs log-like., 10s, agno")
subplot(2,4,2)
PyPlot.plot(1:n,rescale(γ100l1),"o",color="C0")
xlabel("idx")
title("rs ampl., 10s, agno")

subplot(2,4,3)
PyPlot.plot(ks./(T2*τ),rescale(L230l1),color="C1")
xlabel("f")
title("rs log-like., 23s, agno")
subplot(2,4,4)
PyPlot.plot(1:n,rescale(γ230l1),"o",color="C1")
xlabel("idx")
title("rs ampl., 23s, agno")

subplot(2,4,5)
PyPlot.plot(ks./(T1*τ),rescale(L100l1_Aa),color="C2")
xlabel("f")
title("rs log-like., 10s, non-agno")
subplot(2,4,6)
PyPlot.plot(1:n,rescale(γ100l1_Aa),"o",color="C2")
xlabel("idx")
title("rs ampl., 10s, non-agno")

subplot(2,4,7)
PyPlot.plot(ks./(T2*τ),rescale(L230l1_Aa),color="C3")
xlabel("f")
title("rs log-like., 23s, non-agno")
subplot(2,4,8)
PyPlot.plot(1:n,rescale(γ230l1_Aa),"o",color="C3")
xlabel("idx")
title("rs ampl., 23s, non-agno")




L100l0 = readdlm("test-nonagnostic/l0/Agnostic_l0_100_freq",'\t')
L230l0 = readdlm("test-nonagnostic/l0/Agnostic_l0_230_freq",'\t')
L100l0_Aa = readdlm("test-nonagnostic/l0/l0_100_freq",'\t')
L230l0_Aa = readdlm("test-nonagnostic/l0/l0_230_freq",'\t')


figure("SALO")

subplot(2,2,1)
m,ik = findmax(rescale(L100l0))
i = ik[1]
k = ik[2]
idx = [1:i-1;i+1:n]
Lmax = vec(maximum(rescale(L100l0)[idx,:],dims=1))
Lmin = vec(minimum(rescale(L100l0)[idx,:],dims=1))
PyPlot.fill([ks;ks[end:-1:1]]./(T1*τ),[Lmax;Lmin[end:-1:1]],color="gray")
PyPlot.plot(ks./(T1*τ),rescale(L100l0)[i,:],color="C0",label="l = $i")
xlabel("f")
title("rs log-like., 10s, agno")
legend()

subplot(2,2,2)
m,ik = findmax(rescale(L230l0))
i = ik[1]
k = ik[2]
idx = [1:i-1;i+1:n]
Lmax = vec(maximum(rescale(L230l0)[idx,:],dims=1))
Lmin = vec(minimum(rescale(L230l0)[idx,:],dims=1))
PyPlot.fill([ks;ks[end:-1:1]]./(T2*τ),[Lmax;Lmin[end:-1:1]],color="gray")
PyPlot.plot(ks./(T2*τ),rescale(L230l0)[i,:],color="C1",label="l = $i")
xlabel("f")
title("rs log-like., 23s, agno")
legend()

subplot(2,2,3)
m,ik = findmax(rescale(L100l0_Aa))
i = ik[1]
k = ik[2]
idx = [1:i-1;i+1:n]
Lmax = vec(maximum(rescale(L100l0_Aa)[idx,:],dims=1))
Lmin = vec(minimum(rescale(L100l0_Aa)[idx,:],dims=1))
PyPlot.fill([ks;ks[end:-1:1]]./(T1*τ),[Lmax;Lmin[end:-1:1]],color="gray")
PyPlot.plot(ks./(T1*τ),rescale(L100l0_Aa)[i,:],color="C2",label="l = $i")
xlabel("f")
title("rs log-like., 10s, non-agno")
legend()

subplot(2,2,4)
m,ik = findmax(rescale(L230l0_Aa))
i = ik[1]
k = ik[2]
idx = [1:i-1;i+1:n]
Lmax = vec(maximum(rescale(L230l0_Aa)[idx,:],dims=1))
Lmin = vec(minimum(rescale(L230l0_Aa)[idx,:],dims=1))
PyPlot.fill([ks;ks[end:-1:1]]./(T2*τ),[Lmax;Lmin[end:-1:1]],color="gray")
PyPlot.plot(ks./(T2*τ),rescale(L230l0_Aa)[i,:],color="C3",label="l = $i")
xlabel("f")
title("rs log-like., 23s, non-agno")
legend()

