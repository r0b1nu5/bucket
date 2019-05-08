using DelimitedFiles,PyPlot,Statistics,Distributions

include("reservoir.jl")

T_train = 1000
T_simu = 200 

# Building reservoir
N = 1000
m = 10*N
rho = 1.
Adj = A_gen(N,m,rho)
sig = .01
Win = Win_gen(1,N,sig)
a = 1.
xi = 1.
beta = 1.
dT = 1

##=
# Generating unknown system
x0 = rand(2,1)
A = [0 1.;-.1 0]
C = [0 1.]
sig = 0. # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

tt = 10. .+ .1*(1:T_train+dT)
X = Array{Float64,1}()
for t in tt
	push!(X,(C*exp(A*t)*x0)[1]+rand(Normal(0,sig)))
end
ut = [transpose(X[1:end-dT]);]
st = [transpose(X[dT+1:end]);]

Wout,c,rs = reservoir_training((ut,st),Adj,Win,a,xi,beta)

r0 = rs[:,end]
X = Array{Float64,1}()
for t in 10. + .1*T_train .+ .1*(1:dT+T_simu)
	push!(X,(C*exp(A*t)*x0)[1]+rand(Normal(0,sig)))
end
us = [transpose(X[1:T_simu]);]
sc = X[dT+1:end]

ss,rs = reservoir_prediction(us,Wout,c,r0,Adj,Win,a,xi)

PyPlot.plot(sc,"--")
PyPlot.plot(vec(ss),"o")
xlabel("t")
ylabel("y")

# =#
#=
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

# =#




 #= 
# ================= DOES NOT WORK IF INPUT IS TIME =======================


# Generating unknown system
x0 = rand(2,1)
A = [0 1.;-.1 0]
C = [0 1.]
sig = 0. # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

tt = 10*rand(1,T_train)
mtt = mean(tt)
vtt = sqrt(mean((tt .- mtt).^2))
ut = (tt .- mtt)./vtt

Yt = Array{Float64,2}(undef,0,2)
for t in tt
	global A,C,Yt
	Yt = [Yt;C*exp(A*t)]
end
yt = Array{Float64,2}(transpose(Yt*x0 + rand(Normal(0.,sig),T_train)))
myt = sum(yt,dims=2)./T_train
vyt = sqrt.(sum((yt - repeat(myt,1,T_train)).^2,dims=2))
yt = (yt - repeat(myt,1,T_train))./repeat(vyt,1,T_train)

Wout,c,rs = reservoir_training((tt,yt),Adj,Win,a,xi,beta)

r0 = rs[:,end]
ts = 10 .+ 20*rand(1,T_simu)

ys,rs = reservoir_prediction(ts,Wout,c,r0,Adj,Win,a,xi)

yc = Array{Float64,1}()
for t in LinRange(10,30,1000)
	global yc,C,A,x0
	push!(yc,(C*exp(A*t)*x0)[1])
end
yc = (yc .- myt)./vyt

PyPlot.plot(LinRange(10,30,1000),yc,"--")
PyPlot.plot(vec(ts),vec(ys),"o")
xlabel("t")
ylabel("y")

# =#


