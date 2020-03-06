include("reservoir.jl")
include("lorentz.jl")

Ts = [1000,2000,5000,10000,20000,50000,100000,200000,500000]
Ts = [1000000,]

 #=
xs1 = lorentz(rand(3),maximum(Ts)+1000)[:,1001:end]
xs2 = lorentz(rand(3),maximum(Ts)+1000)[:,1001:end]
xs3 = lorentz(rand(3),maximum(Ts)+1000)[:,1001:end]
# =#
# #= 
xs1 = readdlm("data1/xs1.csv",',')
xs2 = readdlm("data1/xs2.csv",',')
xs3 = readdlm("data1/xs3.csv",',')
# =#

N = 500
m = 10*N
rho = 1.
A = A_gen(N,m,rho)
dT = 10

sig = 1.
Win = Win_gen(1,N,sig)

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

# #=
for T in Ts
	W1,c1,r1 = reservoir_training((xs1[[1,],1:T],xs1[2:3,1:T]),A,Win,a,xi,dT,beta)
	W2,c2,r2 = reservoir_training((xs2[[1,],1:T],xs2[2:3,1:T]),A,Win,a,xi,dT,beta)
	W3,c3,r3 = reservoir_training((xs3[[1,],1:T],xs3[2:3,1:T]),A,Win,a,xi,dT,beta)

	push!(Wouts1,W1)
	push!(couts1,c1)
	push!(Wouts2,W2)
	push!(couts2,c2)
	push!(Wouts3,W3)
	push!(couts3,c3)

	push!(noW12,norm(W1 - W2))
	push!(noW23,norm(W2 - W3))
	push!(noW13,norm(W1 - W3))
end

# #=
PyPlot.plot(Ts,noW12,"o")
PyPlot.plot(Ts,noW23,"o")
PyPlot.plot(Ts,noW13,"o")
# =#
# #=
writedlm("data1/Ts.csv",Ts,',')
writedlm("data1/noW12.csv",noW12,',')
writedlm("data1/noW23.csv",noW23,',')
writedlm("data1/noW13.csv",noW13,',')
# =#

# =#
