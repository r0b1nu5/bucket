using DelimitedFiles,PyPlot

include("kernel.jl")

T_train = 5000
T_simu = 50
res = 1
degree = 2
rho = 2.
#gamma = 10000.
gamma = .001

x = readdlm("test/lorentz_100000.csv",',')
	
xt = x[:,res*(1:T_train)]
yt = x[:,res*(2:T_train+1)]-x[:,res*(1:T_train)]

#c,phit = polynomial_kernel_training((xt,yt),degree,gamma)
c = gaussian_kernel_training((xt,yt),rho,gamma)

xs = x[:,res*((T_train-T_simu):(T_train+T_simu))]
x0 = x[:,[res*T_train,]]
ys = Array{Float64,2}(undef,3,0)
# #=
for t in 1:T_simu
#	x0 = polynomial_kernel_prediction(x0,phit,c,degree)
	x0 = gaussian_kernel_prediction(x0,xt,c,rho)
	global ys,x0
	ys = [ys x0]
end
# =#
#ys = polynomial_kernel_prediction(xs,phit,c,degree)
#ys = gaussian_kernel_prediction(xs,xt,c,rho)

#yc = x[:,res*((T_train-T_simu+1):(T_train+T_simu+1))]-x[:,res*((T_train-T_simu):(T_train+T_simu))]
yc = x[:,res*(T_train:T_train+T_simu)]-x[:,res*(T_train+1:T_train+T_simu+1)]

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




