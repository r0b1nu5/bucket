using PyPlot, DelimitedFiles

cols = [(243,111,33)./255,(0,145,194)./255,(161,0,202)./255]
gra = (.8,.8,.8,1.)

ntw = "ntw20"
xy = readdlm("data_melvyn/ntw20/ntw20_xy.csv",',')
A = readdlm("data_melvyn/ntw20/ntw20_A.csv",',')

n = 20
ks = 1:40
τ = .1
T1 = 100
T2 = 230
T3 = 390

function rescale(x::Union{Matrix{Float64}, Vector{Float64}})
	return (abs.(x) .- minimum(abs.(x)))./(maximum(abs.(x)) - minimum(abs.(x)))
end

L100l1 = vec(readdlm("test-nonagnostic/l1/Agnostic_l1_100_freq",','))
γ100l1 = vec(readdlm("test-nonagnostic/l1/Agnostic_l1_100_gamma",','))
L230l1 = vec(readdlm("test-nonagnostic/l1/Agnostic_l1_230_freq",','))
γ230l1 = vec(readdlm("test-nonagnostic/l1/Agnostic_l1_230_gamma",','))
Lxxxl1 = vec(readdlm("test-nonagnostic/l1/Agnostic_l1_$(T3)_freq",','))
γxxxl1 = vec(readdlm("test-nonagnostic/l1/Agnostic_l1_$(T3)_gamma",','))

L100l1_Aa = vec(readdlm("test-nonagnostic/l1/l1_100_freq",','))
γ100l1_Aa = vec(readdlm("test-nonagnostic/l1/l1_100_gamma",','))
L230l1_Aa = vec(readdlm("test-nonagnostic/l1/l1_230_freq",','))
γ230l1_Aa = vec(readdlm("test-nonagnostic/l1/l1_230_gamma",','))
Lxxxl1_Aa = vec(readdlm("test-nonagnostic/l1/l1_$(T3)_freq",','))
γxxxl1_Aa = vec(readdlm("test-nonagnostic/l1/l1_$(T3)_gamma",','))

figure("SALO-relax")

#subplot(1,4,1)
subplot(2,2,1)
for i in 1:n-1
	for j in i+1:n
		if A[i,j] > 1e-2
			PyPlot.plot(xy[[i,j],1],-xy[[i,j],2],"k")
		end
	end
end
PyPlot.plot(xy[:,1],-xy[:,2],"ok",markersize=6.)
PyPlot.plot(xy[20,1],-xy[20,2],"o",color=cols[1],markersize=8.)
PyPlot.plot(xy[14,1],-xy[14,2],"o",color=cols[2],markersize=8.)
PyPlot.plot(xy[4,1],-xy[4,2],"o",color=cols[3],markersize=8.)
axis([2*minimum(xy[:,1]),2*maximum(xy[:,1]),2*minimum(xy[:,2]),2*maximum(xy[:,2])])
xticks([])
yticks([])

#subplot(2,8,3)
subplot(4,4,3)
PyPlot.plot([.08,.08],[-.1,1.1],"--",color="C7")
PyPlot.plot([0;ks]./(T1*τ),[0;rescale(L100l1)],color=cols[2])
axis([0/(T1*τ),40/(T1*τ),-.1,1.1])
xlabel("f")
ylabel("rescaled \nlog-likelihood")
title("10s, A unknown")
#subplot(2,8,11)
subplot(4,4,7)
PyPlot.plot(1:n,rescale(γ100l1),"o",color=cols[2])
axis([0,n+1,-.1,1.1])
xlabel("idx")
ylabel("rescaled \namplitude")
#title("10s, A unknown")

#subplot(2,8,5)
subplot(4,4,9)
PyPlot.plot([.08,.08],[-.1,1.1],"--",color="C7")
PyPlot.plot(ks./(T2*τ),rescale(L230l1),color=cols[3])
axis([1/(T2*τ),40/(T2*τ),-.1,1.1])
xlabel("f")
ylabel("rescaled \nlog-likelihood")
title("23s, A unknown")
#subplot(2,8,13)
subplot(4,4,13)
PyPlot.plot(1:n,rescale(γ230l1),"o",color=cols[3])
axis([0,n+1,-.1,1.1])
xlabel("idx")
ylabel("rescaled \namplitude")
#title("23s, A unknown")

#subplot(2,8,4)
subplot(4,4,4)
PyPlot.plot([.08,.08],[-.1,1.1],"--",color="C7")
PyPlot.plot([0;ks]./(T1*τ),[0;rescale(L100l1_Aa)],color=cols[1])
axis([0/(T1*τ),40/(T1*τ),-.1,1.1])
xlabel("f")
title("10s, A known")
#subplot(2,8,12)
subplot(4,4,8)
PyPlot.plot(1:n,rescale(γ100l1_Aa),"o",color=cols[1])
axis([0,n+1,-.1,1.1])
xlabel("idx")
#title("10s, A known")

#subplot(2,8,6)
subplot(4,4,10)
PyPlot.plot([.08,.08],[-.1,1.1],"--",color="C7")
PyPlot.plot(ks./(T2*τ),rescale(L230l1_Aa),color=cols[1])
axis([1/(T2*τ),40/(T2*τ),-.1,1.1])
xlabel("f")
title("23s, A known")
#subplot(2,8,14)
subplot(4,4,14)
PyPlot.plot(1:n,rescale(γ230l1_Aa),"o",color=cols[1])
axis([0,n+1,-.1,1.1])
xlabel("idx")
#title("23s, A known")

#subplot(2,8,7)
subplot(4,4,12)
PyPlot.plot([.08,.08],[-.1,1.1],"--",color="C7")
PyPlot.plot(ks./(T3*τ),rescale(Lxxxl1_Aa),color=cols[1])
axis([1/(T3*τ),40/(T3*τ),-.1,1.1])
xlabel("f")
title("$(Int64(T3*τ))s, A known")
#subplot(2,8,15)
subplot(4,4,16)
PyPlot.plot(1:n,rescale(γxxxl1_Aa),"o",color=cols[1])
axis([0,n+1,-.1,1.1])
xlabel("idx")
#title("$(Int64(T3*τ))s, A known")

#subplot(2,8,8)
subplot(4,4,11)
PyPlot.plot([.08,.08],[-.1,1.1],"--",color="C7")
PyPlot.plot(ks./(T3*τ),rescale(Lxxxl1),color=cols[1])
axis([1/(T3*τ),40/(T3*τ),-.1,1.1])
xlabel("f")
title("$(Int64(T3*τ))s, A unknown")
#subplot(2,8,16)
subplot(4,4,15)
PyPlot.plot(1:n,rescale(γxxxl1),"o",color=cols[1])
axis([0,n+1,-.1,1.1])
xlabel("idx")
#title("$(Int64(T3*τ))s, A unknown")




L100l0 = readdlm("test-nonagnostic/l0/Agnostic_l0_100_freq",'\t')
L100l0_Aa = readdlm("test-nonagnostic/l0/l0_100_freq",'\t')
L230l0 = readdlm("test-nonagnostic/l0/Agnostic_l0_230_freq",'\t')
L230l0_Aa = readdlm("test-nonagnostic/l0/l0_230_freq",'\t')
Lxxxl0 = readdlm("test-nonagnostic/l0/Agnostic_l0_$(T3)_freq",'\t')
Lxxxl0_Aa = readdlm("test-nonagnostic/l0/l0_$(T3)_freq",'\t')


figure("SALO")

#subplot(1,4,1)
subplot(2,2,1)
for i in 1:n-1
	for j in i+1:n
		if A[i,j] > 1e-2
			PyPlot.plot(xy[[i,j],1],-xy[[i,j],2],"k")
		end
	end
end
PyPlot.plot(xy[:,1],-xy[:,2],"ok",markersize=6.)
PyPlot.plot(xy[20,1],-xy[20,2],"o",color=cols[1],markersize=8.)
PyPlot.plot(xy[8,1],-xy[8,2],"o",color=cols[2],markersize=8.)
PyPlot.plot(xy[19,1],-xy[19,2],"o",color=cols[3],markersize=8.)
axis([2*minimum(xy[:,1]),2*maximum(xy[:,1]),2*minimum(xy[:,2]),2*maximum(xy[:,2])])
xticks([])
yticks([])

#subplot(2,4,2)
subplot(4,2,2)
m,ik = findmax(rescale(L100l0))
i = ik[1]
k = ik[2]
idx = [1:i-1;i+1:n]
Lmax = vec(maximum(rescale(L100l0)[idx,:],dims=1))
Lmin = vec(minimum(rescale(L100l0)[idx,:],dims=1))
PyPlot.fill([ks;ks[end:-1:1]]./(T1*τ),[Lmax;Lmin[end:-1:1]],color=gra)
PyPlot.plot([.08,.08],[-.1,1.1],"--",color="C7")
PyPlot.plot(ks./(T1*τ),rescale(L100l0)[i,:],color=cols[2],label="l = $i")
axis([1/(T1*τ),40/(T1*τ),-.1,1.1])
xlabel("f")
ylabel("rescaled \nlog-likelihood")
title("10s, A unknown")
legend()

#subplot(2,4,3)
subplot(4,2,5)
m,ik = findmax(rescale(L230l0))
i = ik[1]
k = ik[2]
idx = [1:i-1;i+1:n]
Lmax = vec(maximum(rescale(L230l0)[idx,:],dims=1))
Lmin = vec(minimum(rescale(L230l0)[idx,:],dims=1))
PyPlot.fill([ks;ks[end:-1:1]]./(T2*τ),[Lmax;Lmin[end:-1:1]],color=gra)
PyPlot.plot([.08,.08],[-.1,1.1],"--",color="C7")
PyPlot.plot(ks./(T2*τ),rescale(L230l0)[i,:],color=cols[3],label="l = $i")
axis([1/(T2*τ),40/(T2*τ),-.1,1.1])
xlabel("f")
ylabel("rescaled \nlog-likelihood")
title("23s, A unknown")
legend()

#subplot(2,4,6)
subplot(4,2,4)
m,ik = findmax(rescale(L100l0_Aa))
i = ik[1]
k = ik[2]
idx = [1:i-1;i+1:n]
Lmax = vec(maximum(rescale(L100l0_Aa)[idx,:],dims=1))
Lmin = vec(minimum(rescale(L100l0_Aa)[idx,:],dims=1))
PyPlot.fill([ks;ks[end:-1:1]]./(T1*τ),[Lmax;Lmin[end:-1:1]],color=gra)
PyPlot.plot([.08,.08],[-.1,1.1],"--",color="C7")
PyPlot.plot(ks./(T1*τ),rescale(L100l0_Aa)[i,:],color=cols[1],label="l = $i")
axis([1/(T1*τ),40/(T1*τ),-.1,1.1])
xlabel("f")
ylabel("rescaled \nlog-likelihood")
title("10s, A known")
legend()

#subplot(2,4,7)
subplot(4,2,7)
m,ik = findmax(rescale(L230l0_Aa))
i = ik[1]
k = ik[2]
idx = [1:i-1;i+1:n]
Lmax = vec(maximum(rescale(L230l0_Aa)[idx,:],dims=1))
Lmin = vec(minimum(rescale(L230l0_Aa)[idx,:],dims=1))
PyPlot.fill([ks;ks[end:-1:1]]./(T2*τ),[Lmax;Lmin[end:-1:1]],color=gra)
PyPlot.plot([.08,.08],[-.1,1.1],"--",color="C7")
PyPlot.plot(ks./(T2*τ),rescale(L230l0_Aa)[i,:],color=cols[1],label="l = $i")
axis([1/(T2*τ),40/(T2*τ),-.1,1.1])
xlabel("f")
ylabel("rescaled \nlog-likelihood")
title("23s, A known")
legend()

#subplot(2,4,4)
subplot(4,2,6)
m,ik = findmax(rescale(Lxxxl0))
i = ik[1]
k = ik[2]
idx = [1:i-1;i+1:n]
Lmax = vec(maximum(rescale(Lxxxl0)[idx,:],dims=1))
Lmin = vec(minimum(rescale(Lxxxl0)[idx,:],dims=1))
PyPlot.fill([ks;ks[end:-1:1]]./(T3*τ),[Lmax;Lmin[end:-1:1]],color=gra)
PyPlot.plot([.08,.08],[-.1,1.1],"--",color="C7")
PyPlot.plot(ks./(T3*τ),rescale(Lxxxl0)[i,:],color=cols[1],label="l = $i")
axis([1/(T3*τ),40/(T3*τ),-.1,1.1])
xlabel("f")
title("$(Int64(T3*τ))s, A unknown")
legend()

#subplot(2,4,8)
subplot(4,2,8)
m,ik = findmax(rescale(Lxxxl0_Aa))
i = ik[1]
k = ik[2]
idx = [1:i-1;i+1:n]
Lmax = vec(maximum(rescale(Lxxxl0_Aa)[idx,:],dims=1))
Lmin = vec(minimum(rescale(Lxxxl0_Aa)[idx,:],dims=1))
PyPlot.fill([ks;ks[end:-1:1]]./(T3*τ),[Lmax;Lmin[end:-1:1]],color=gra)
PyPlot.plot([.08,.08],[-.1,1.1],"--",color="C7")
PyPlot.plot(ks./(T3*τ),rescale(Lxxxl0_Aa)[i,:],color=cols[1],label="l = $i")
axis([1/(T3*τ),40/(T3*τ),-.1,1.1])
xlabel("f")
title("$(Int64(T3*τ))s, A known")
legend()
