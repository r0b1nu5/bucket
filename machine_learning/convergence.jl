include("reservoir.jl")
include("lorenz.jl")

Ts = [100,200,500,1000,2000,5000,10000,20000,50000,100000]
Ts = Array(5001:500:10001)

 #=
xs1 = lorenz(rand(3),maximum(Ts)+1000)[:,1001:end]
xs2 = lorenz(rand(3),maximum(Ts)+1000)[:,1001:end]
xs3 = lorenz(rand(3),maximum(Ts)+1000)[:,1001:end]
# =#
# #= 
xs1 = readdlm("data1/xs1.csv",',')
xs2 = readdlm("data1/xs2.csv",',')
xs3 = readdlm("data1/xs3.csv",',')
# =#

n = size(xs1)[1]

N = 500
m = 10*N
rho = 1.
A = A_gen(N,m,rho)
dT = 1

sig = 1.
Win = Win_gen(n,N,sig)

a = 1.
xi = 1.
beta = .01

Wouts1 = Array{Array{Float64,2},1}()
Wouts2 = Array{Array{Float64,2},1}()
Wouts3 = Array{Array{Float64,2},1}()
couts1 = Array{Array{Float64,1},1}()
couts2 = Array{Array{Float64,1},1}()
couts3 = Array{Array{Float64,1},1}()
noW12 = Array{Float64,1}()
noW23 = Array{Float64,1}()
noW13 = Array{Float64,1}()
cauchy1 = Array{Float64,1}()
cauchy2 = Array{Float64,1}()
cauchy3 = Array{Float64,1}()
svd1 = Array{Float64,2}(undef,minimum(size(xs1)),0)
svd2 = Array{Float64,2}(undef,minimum(size(xs1)),0)
svd3 = Array{Float64,2}(undef,minimum(size(xs1)),0)

# #=
for T in Ts
	W1,c1,r1 = reservoir_training((xs1[:,1:T-1],xs1[:,2:T]),A,Win,a,xi,dT,beta)
	W2,c2,r2 = reservoir_training((xs2[:,1:T-1],xs2[:,2:T]),A,Win,a,xi,dT,beta)
	W3,c3,r3 = reservoir_training((xs3[:,1:T-1],xs3[:,2:T]),A,Win,a,xi,dT,beta)

	push!(Wouts1,W1)
	push!(couts1,c1)
	push!(Wouts2,W2)
	push!(couts2,c2)
	push!(Wouts3,W3)
	push!(couts3,c3)

	push!(noW12,norm(W1 - W2))
	push!(noW23,norm(W2 - W3))
	push!(noW13,norm(W1 - W3))

	if T > minimum(Ts)
		push!(cauchy1,norm(W1-Wouts1[end-1]))
		push!(cauchy2,norm(W2-Wouts2[end-1]))
		push!(cauchy3,norm(W3-Wouts3[end-1]))
	end

	svd1 = [svd1 svd(W1).S]
	svd2 = [svd2 svd(W2).S]
	svd3 = [svd3 svd(W3).S]

	global svd1,svd2,svd3
end

# #=
figure()
PyPlot.plot(Ts,noW12,"-o",label="||W1-W2||")
PyPlot.plot(Ts,noW23,"-o",label="||W2-W3||")
PyPlot.plot(Ts,noW13,"-o",label="||W1-W3||")
legend()
xlabel("Training time")
ylabel("Norm")

figure()
PyPlot.plot(Ts[2:end],cauchy1,"-o",label="||W1(T)-W1(T-δT)||")
PyPlot.plot(Ts[2:end],cauchy2,"-o",label="||W2(T)-W2(T-δT)||")
PyPlot.plot(Ts[2:end],cauchy3,"-o",label="||W3(T)-W3(T-δT)||")
legend()
xlabel("Training time")
ylabel("Norm")

figure()
subplot(3,1,1)
PyPlot.plot(Ts,svd1[1,:],"-o")
PyPlot.plot(Ts,svd2[1,:],"-o")
PyPlot.plot(Ts,svd3[1,:],"-o")
ylabel("1st SV")
subplot(3,1,2)
PyPlot.plot(Ts,svd1[2,:],"-o")
PyPlot.plot(Ts,svd2[2,:],"-o")
PyPlot.plot(Ts,svd3[2,:],"-o")
ylabel("2nd SV")
subplot(3,1,3)
PyPlot.plot(Ts,svd1[3,:],"-o")
PyPlot.plot(Ts,svd2[3,:],"-o")
PyPlot.plot(Ts,svd3[3,:],"-o")
ylabel("3rd SV")

# =#
# #=
writedlm("data1/Ts.csv",Ts,',')
writedlm("data1/noW12.csv",noW12,',')
writedlm("data1/noW23.csv",noW23,',')
writedlm("data1/noW13.csv",noW13,',')
# =#

# =#
