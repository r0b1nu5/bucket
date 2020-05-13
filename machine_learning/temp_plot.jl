using PyPlot, DelimitedFiles

ids = ["_pj1","_pj2","_pj3"]
thrs = [.1,.05,.02]
ls = [1.50255,1.37446,0.90566]
data = ["(16.0,45.92,4.0)", "(16.0,40.0,4.0)", "(10.0,28.0,2.66)"]

N = 100
Tbs = zeros(3,3,N)

for i in 1:3
	for j in 1:3
		for k in 1:N
			global Tbs,thrs
			Tbs[i,j,k] = readdlm("data1/Tb"*ids[i]*"_$(k).1.1_vs_N_vs_Tt_thr$(thrs[j]).csv",',')[1]
		end
		PyPlot.plot(1/ls[i]+(j-2)*2e-3,quantile(Tbs[i,j,:],.5),"-o",color="C$(j-1)")
		PyPlot.plot((1/ls[i]+(j-2)*2e-3)*[1,1],[quantile(Tbs[i,j,:],.25),quantile(Tbs[i,j,:],.75)],"-",color="C$(j-1)")
#		PyPlot.plot((1/ls[i]+(j-2)*2e-3)*[1,1],[quantile(Tbs[i,j,:],0.),quantile(Tbs[i,j,:],1.)],"x",color="C$(j-1)")
	end
end

xlabel("1/λ")
ylabel("Breaktime")
title("Lorenz: N = 2000, Tt = 2001")



