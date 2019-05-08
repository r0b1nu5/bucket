using DelimitedFiles,PyPlot,Statistics

include("kernel.jl")

#kernel_type = "polynomial"
kernel_type = "gaussian"

T_train = 5000
T_simu = 200
res = 50
degree = 2
rho = 1.
gamma = 1e-9
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

# #=
figure()
PyPlot.plot(yc[1,:])
PyPlot.plot(ys[1,:],"--")
xlabel("t")
ylabel("x")
figure()
PyPlot.plot(yc[2,:])
PyPlot.plot(ys[2,:],"--")
xlabel("t")
ylabel("y")
figure()
PyPlot.plot(yc[3,:])
PyPlot.plot(ys[3,:],"--")
xlabel("t")
ylabel("z")
figure()
PyPlot.semilogy(vec(sum((yc-ys).^2,dims=1)))
xlabel("t")
ylabel("error")
# =#




