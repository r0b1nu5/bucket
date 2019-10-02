using PyPlot, LinearAlgebra, Statistics, FFTW

figure(69)

subplot(222)
for e in 1:9
	PyPlot.semilogy([-.0001,.015],[(10.)^e,(10.)^e],"--k",linewidth=1)
end

subplot(224)
for e in -1:5
	PyPlot.semilogy([-.0001,.015],[(10.)^e,(10.)^e],"--k",linewidth=1)
end

Xs = readdlm("data2/backup2_ntw_forced_0.009_100000_0.1.csv",',')
X = Xs[6,:] .- mean(Xs[6,:])
co = "C0"
sty = "-"

subplot(221)
PyPlot.plot(dt*(1:T),X,sty,color=co)

fX = fft(X)
nfX = norm.(fX)
fs = (0:T-1)./(dt*T)
subplot(222)
PyPlot.semilogy(fs[2:end],nfX[2:end],sty,color=co)

R = Array{Float64,1}()
for i in 1:Int(T/10)
	push!(R,mean(X[1:end-i].*X[1+i:end]))
end
subplot(223)
PyPlot.plot((1:Int(T/10))*dt,R,sty,color=co)

Df = 10
mfX = Array{Float64,1}()
for j in 1:T
	push!(mfX,nfX[j]/median([nfX[max(j-Df,1):j-2];nfX[j+2:min(j+Df,T)]]))
end
subplot(224)
PyPlot.semilogy(fs[2:end],mfX[2:end],sty,color=co)


Xs = readdlm("data2/backup1_ntw_forced_0.009_100000_0.1.csv",',')
X = Xs[6,:] .- mean(Xs[6,:])
co = "C1"
sty = "-"

subplot(221)
PyPlot.plot(dt*(1:T),X .+ 4,sty,color=co)

fX = fft(X)
nfX = norm.(fX)
fs = (0:T-1)./(dt*T)
subplot(222)
PyPlot.semilogy(fs[2:end],nfX[2:end] * 100,sty,color=co)

R = Array{Float64,1}()
for i in 1:Int(T/10)
	push!(R,mean(X[1:end-i].*X[1+i:end]))
end
subplot(223)
PyPlot.plot((1:Int(T/10))*dt,R .- 1.3,sty,color=co)

Df = 10
mfX = Array{Float64,1}()
for j in 1:T
	push!(mfX,nfX[j]/median([nfX[max(j-Df,1):j-2];nfX[j+2:min(j+Df,T)]]))
end
subplot(224)
PyPlot.semilogy(fs[2:end],mfX[2:end] * 100,sty,color=co)


Xs = readdlm("data2/backup3_ntw_forced_0.009_100000_0.1.csv",',')
X = Xs[6,:] .- mean(Xs[6,:])
X ./= 5
co = "C2"
sty = "-"

subplot(221)
PyPlot.plot(dt*(1:T),X .+ 10,sty,color=co)

fX = fft(X)
nfX = norm.(fX)
fs = (0:T-1)./(dt*T)
subplot(222)
PyPlot.semilogy(fs[2:end],nfX[2:end] * 10000,sty,color=co)

R = Array{Float64,1}()
for i in 1:Int(T/10)
	push!(R,mean(X[1:end-i].*X[1+i:end]))
end
subplot(223)
PyPlot.plot((1:Int(T/10))*dt,R .+ .5,sty,color=co)

Df = 10
mfX = Array{Float64,1}()
for j in 1:T
	push!(mfX,nfX[j]/median([nfX[max(j-Df,1):j-2];nfX[j+2:min(j+Df,T)]]))
end
subplot(224)
PyPlot.semilogy(fs[2:end],mfX[2:end] * 10000,sty,color=co)


subplot(222)
axis([-.0001,.015,1e1,1e9])

subplot(224)
axis([-.0001,.015,.1,2e5])



