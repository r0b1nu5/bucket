using DelimitedFiles,PyPlot,Statistics

include("reservoir.jl")

T_train = 30
T_simu = 200 
res = 5

#x = readdlm("test/lorentz_100000.csv",',')

n = size(x)[1]
N = 500
m = 10*N
rho = 1.
A = A_gen(N,m,rho)
sig = 1.
Win = Win_gen(1,N,sig)
a = 1.
xi = 1.
beta = .01
T0 = 200

ut = x[[1,],T0.+res*(1:T_train)]
mut = mean(ut)
vut = sqrt(mean((ut .- mut).^2))
ut = (ut .- mut)./vut

st = x[2:3,T0.+res*(1:T_train)]
mst = sum(st,dims=2)./T_train
vst = sqrt.(sum((st - repeat(mst,1,T_train)).^2,dims=2))
st = (st - repeat(mst,1,T_train))./repeat(vst,1,T_train)

Wout,c,rs = reservoir_training((ut,st),A,Win,a,xi,beta)

r0 = rs[:,end]
us = x[[1,],T0.+res*((T_train+1):(T_train+T_simu))]
us = (us .- mut)./vut

ss,rs = reservoir_prediction(us,Wout,c,r0,A,Win,a,xi)

ys = [us;ss]
yc = x[:,T0.+res*((T_train+1):(T_train+T_simu))]
yc = (yc - repeat([mut;mst],1,T_simu))./repeat([vut;vst],1,T_simu)

figure()
PyPlot.plot(yc[1,:])
PyPlot.plot(ys[1,:],"--")
figure()
PyPlot.plot(yc[2,:])
PyPlot.plot(ys[2,:],"--")
figure()
PyPlot.plot(yc[3,:])
PyPlot.plot(ys[3,:],"--")
figure()
PyPlot.semilogy(vec(sum((yc-ys).^2,dims=1)))




