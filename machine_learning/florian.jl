using PyPlot,Statistics,Distributions

N = 100
l = .001
t0 = 10.
T = 200

A = [0 1;-.1 0]
C = [0 1.]

x0 = rand(2)

ts = 10*rand(N)
xs = Array{Float64,2}(undef,2,0)
ys = Array{Float64,1}()
for t in ts
	xs = [xs exp(A*t)*x0]
	push!(ys,(C*xs[:,end])[1] + rand(Normal(0.,.1)))
	global xs
end

K = zeros(N,N)
for i in 1:N-1
	K[i,i] = (C*exp(A*ts[i])*transpose(exp(A*ts[i]))*transpose(C))[1]
	for j in i+1:N
		K[i,j] = (C*exp(A*ts[i])*transpose(exp(A*ts[j]))*transpose(C))[1]
		K[j,i] = K[i,j]
	end
	global K
end

Opt = inv(l*diagm(0 => ones(N)) + K)*ys

yz = Array{Float64,1}()
yr = Array{Float64,1}()

for t in t0+1:t0+T
	k = Array{Float64,1}()
	for i in 1:N
		push!(k,(C*exp(A*t)*transpose(exp(A*ts[i]))*transpose(C))[1])
	end
	push!(yz,transpose(k)*Opt)
	push!(yr,(C*exp(A*t)*x0)[1])
	global x0,yr,yz
end

PyPlot.plot(yr,"--")
PyPlot.plot(yz,"o")


