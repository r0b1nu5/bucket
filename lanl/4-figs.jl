using PyPlot, LinearAlgebra, Statistics, FFTW

# #=
co = "C0"
Xs = generate_forced_time_series("ntw",L,m,d,(a,f,phi),T,dt,sig)
X = Xs[6,:] .- mean(Xs[6,:])
# =#

figure(70)

subplot(221)
PyPlot.plot(dt*(1:T),X,color=co)

fX = fft(X)
nfX = norm.(fX)
fs = (0:T-1)./(dt*T)
subplot(222)
PyPlot.plot(fs,nfX,color=co)

R = Array{Float64,1}()
for i in 1:Int(T/10)
	push!(R,mean(X[1:end-i].*X[1+i:end]))
end
subplot(223)
PyPlot.plot((1:Int(T/10))*dt,R,color=co)

Df = 10
mfX = Array{Float64,1}()
for j in 1:T
	push!(mfX,nfX[j]/median([nfX[max(j-Df,1):j-2];nfX[j+2:min(j+Df,T)]]))
end
subplot(224)
PyPlot.plot(fs,mfX,color=co)



