using DelimitedFiles, FFTW

include("generate_time_series.jl")

#regime = "low"
#regime = "mid"
regime = "high"

n = 5

L = readdlm("data/ntw5_lap_mat.csv",',')
m = vec(readdlm("data/ntw5_m.csv",','))
d = vec(readdlm("data/ntw5_d.csv",','))

a = [.2,0.,0.,0.,0.]
f = zeros(n)
if regime == "low"
	f[1] = .009
	ax1 = [500,1000,-.18,.18]
	ax2 = [0.,.04,1e-1,1e4]
elseif regime == "mid"
	f[1] = .7
	ax1 = [1000,1050,-.18,.18]
	ax2 = [0.,1.,1e-1,1e4]
elseif regime == "high"
	f[1] = 2.
	ax1 = [1000,1020,-.15,.15]
	ax2 = [0.,2.2,1e-2,1e4]
end
phi = [pi/10,0.,0.,0.,0.]

T = 100000
dt = .1

Xs = generate_forced_time_series("ntw5",L,m,d,(a,f,phi),T,dt,.2*m)

figure()

subplot(5,2,1)
PyPlot.plot((1:T).*dt,Xs[6,:],color="C0")
xlabel("t [s]")
ylabel("thd")
axis(ax1)

subplot(5,2,2)
PyPlot.semilogy((0:T-1)./(dt*T),norm.(fft(Xs[6,:])),color="C0")
xlabel("f [2pi/s]")
ylabel("DFT")
axis(ax2)

subplot(5,2,3)
PyPlot.plot((1:T).*dt,Xs[7,:],color="C1")
xlabel("t [s]")
ylabel("thd")
axis(ax1)

subplot(5,2,4)
PyPlot.semilogy((0:T-1)./(dt*T),norm.(fft(Xs[7,:])),color="C1")
xlabel("f [2pi/s]")
ylabel("DFT")
axis(ax2)

subplot(5,2,5)
PyPlot.plot((1:T).*dt,Xs[8,:],color="C2")
xlabel("t [s]")
ylabel("thd")
axis(ax1)

subplot(5,2,6)
PyPlot.semilogy((0:T-1)./(dt*T),norm.(fft(Xs[8,:])),color="C2")
xlabel("f [2pi/s]")
ylabel("DFT")
axis(ax2)

subplot(5,2,7)
PyPlot.plot((1:T).*dt,Xs[9,:],color="C3")
xlabel("t [s]")
ylabel("thd")
axis(ax1)

subplot(5,2,8)
PyPlot.semilogy((0:T-1)./(dt*T),norm.(fft(Xs[9,:])),color="C3")
xlabel("f [2pi/s]")
ylabel("DFT")
axis(ax2)

subplot(5,2,9)
PyPlot.plot((1:T).*dt,Xs[10,:],color="C4")
xlabel("t [s]")
ylabel("thd")
axis(ax1)

subplot(5,2,10)
PyPlot.semilogy((0:T-1)./(dt*T),norm.(fft(Xs[10,:])),color="C4")
xlabel("f [2pi/s]")
ylabel("DFT")
axis(ax2)



