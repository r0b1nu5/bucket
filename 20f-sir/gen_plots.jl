using PyPlot, DelimitedFiles

include("sir.jl")

plot_bound(1)
plot_data("data/I_ref_0229.csv",1)

N = 8500000
start_day = 13
n = Int(readdlm("data/I_ref_0229.csv",',')[start_day+1])

Rref = 2.7

for R0 in [Rref,2.,1.7,1.5,1.35]
	sirs = sir(N,(Rref,R0,0,200,1000000),5.,[(N-n)/N,n/N,0.],start_day,1,1e-6)
end

figure(1)
subplot(2,1,1)
axis([-1,150,-.01,.3])
legend()
xlabel("Days since Feb 29th")
ylabel("Prop. infected [%]")
subplot(2,1,2)
axis([-1,150,-.01,.3])
xlabel("Days since Feb 29th")
ylabel("Prop. infected [%]")

plot_bound(2)
plot_data("data/I_ref_0229.csv",2)
sirs = sir(N,(Rref,Rref,700,700,1000000),5.,[(N-n)/N,n/N,0.],start_day,2,1e-6)

R0 = 1.35
for ti in [1400,1000,700,200]
	sirs = sir(N,(Rref,R0,ti,700,1000000),5.,[(N-n)/N,n/N,0.],start_day,2,1e-6)
end

figure(2)
subplot(2,1,1)
axis([-1,150,-.01,.3])
legend()
xlabel("Days since Feb 29th")
ylabel("Prop. infected [%]")
subplot(2,1,2)
axis([-1,150,-.01,.3])
xlabel("Days since Feb 29th")
ylabel("Prop. infected [%]")

plot_bound(3)
plot_data("data/I_ref_0229.csv",3)
sirs = sir(N,(Rref,Rref,200,700,1000000),5.,[(N-n)/N,n/N,0.],start_day,3,1e-6)

for tt in [2000,1400,1000,700,200]
	sirs = sir(N,(Rref,R0,700,tt,1000000),5.,[(N-n)/N,n/N,0.],start_day,3,1e-6)
end

figure(3)
subplot(2,1,1)
axis([-1,150,-.01,.3])
legend()
xlabel("Days since Feb 29th")
ylabel("Prop. infected [%]")
subplot(2,1,2)
axis([-1,150,-.01,.3])
xlabel("Days since Feb 29th")
ylabel("Prop. infected [%]")





