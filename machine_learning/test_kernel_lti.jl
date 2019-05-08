using DelimitedFiles,PyPlot,Statistics,Distributions

include("kernel.jl")

#kernel_type = "polynomial"
#kernel_type = "gaussian"
kernel_type = "lti"

T = 2
gamma = .001

A = [0 1.;-.1 0]
C = [0 1.]
sig = .01

x0 = rand(2)
tt = 10*rand(T)
Yt = Array{Float64,2}(undef,0,2)
yt = Array{Float64,1}()
for t in tt
	global x0,A,C,Yt,yt,sig
	Yt = [Yt;C*exp(A*t)]
end
yt = Yt*x0 + rand(Normal(0.,sig),T)

K = Yt*transpose(Yt)
Opt = inv(gamma*diagm(0 => ones(T)) + K)*yt

ts = 10 .+ 20*rand(100)
ys = Array{Float64,1}()
for t in ts
	global A,C,ys,tt,Opt
	x = Array{Float64,1}()
	for t2 in tt
		push!(x,(C*exp(A*t)*transpose(exp(A*t2))*transpose(C))[1])
	end
	push!(ys,(transpose(x)*Opt)[1])
end

yc = Array{Float64,1}()
for t in LinRange(10,30,1000)
	global yc,C,A,x0
	push!(yc,(C*exp(A*t)*x0)[1])
end

PyPlot.plot(LinRange(10,30,1000),yc,"--")
PyPlot.plot(ts,ys,"o")


#=


T_train = 5000
T_simu = 200
res = 50
degree = 2
rho = 1.
#gamma = 1e-8
#gamma = 5e-14
T0 = 20000

#x = readdlm("test/lorentz_100000.csv",',')
	
xt = x[[1,],T0.+res*(1:T_train)]
mxt = mean(xt)
vxt = sqrt(mean((xt .- mxt).^2))
xt = (xt .- mxt)./vxt

yt = x[2:3,T0.+res*(1:T_train)]
myt = sum(yt,dims=2)./T_train
vyt = sqrt.(sum((yt - repeat(myt,1,T_train)).^2,dims=2))
yt = (yt - repeat(myt,1,T_train))./repeat(vyt,1,T_train)

xs = Array{Float64,2}(undef,0,0)
ys = Array{Float64,2}(undef,0,0)

if kernel_type == "polynomial"
	global xt,yt,degree,gamma,xs,ys
	
	c,phit = polynomial_kernel_training((xt,yt),degree,gamma)
	
	xs = x[[1,],T0.+res*((T_train+1):(T_train+T_simu))]
	
	ys = polynomial_kernel_prediction(xs,phit,c,degree)
elseif kernel_type == "gaussian"
	c = gaussian_kernel_training((xt,yt),rho,gamma)
	
	xs = x[[1,],T0.+res*((T_train+1):(T_train+T_simu))]
	
	ys = gaussian_kernel_prediction(xs,xt,c,rho)
end

yc = x[:,T0.+res*((T_train+1):(T_train+T_simu))]
ys = [xs;ys]

 #=
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



# =#
