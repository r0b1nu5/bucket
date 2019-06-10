using PyPlot,DelimitedFiles

include("../default_color_cycle.jl")

Es = abs.(readdlm("data/Es_4.csv",','))
vars = readdlm("data/vars_4.csv",',')

figure()
subplot(141)
PyPlot.loglog(Es[1,:],Es[3,:],"x",color=def_col[1])
PyPlot.plot([1e-8,1e2],[1e-8,1e2],"--k")
axis([minimum(Es[[1,3],:]),maximum(Es[[1,3],:]),minimum(Es[[1,3],:]),maximum(Es[[1,3],:])])
xlabel("Predicted expectation (effective data)")
ylabel("Measured expectabtion")

subplot(143)
PyPlot.loglog(vars[1,:],vars[3,:],"x",color=def_col[2])
PyPlot.plot([1e-8,1e2],[1e-8,1e2],"--k")
axis([minimum(vars[[1,3],:]),maximum(vars[[1,3],:]),minimum(vars[[1,3],:]),maximum(vars[[1,3],:])])
xlabel("Predicted variance (effective data)")
ylabel("Measured variance")

ex = 0.
Ex = 0.

vx = 0.
Vx = 0.

for k in [6,]
	global ex,Ex,vx,Vx

	Es = abs.(readdlm("data/Es_$k.csv",','))
	vars = readdlm("data/vars_$k.csv",',')
	
	ex = max(ex,minimum(Es[[1,3],:]))
	Ex = max(Ex,maximum(Es[[1,3],:]))
	vx = max(vx,minimum(vars[[1,3],:]))
	Vx = max(Vx,maximum(vars[[1,3],:]))

	subplot(142)
	PyPlot.loglog(Es[1,:],Es[3,:],"x",color=def_col[1])

	subplot(144)
	PyPlot.loglog(vars[1,:],vars[3,:],"x",color=def_col[2])
end

subplot(142)
PyPlot.plot([1e-8,1e2],[1e-8,1e2],"--k")
axis([ex,Ex,ex,Ex])
xlabel("Predicted expectation (theoretical data)")
ylabel("Measured expectabtion")

subplot(144)
PyPlot.plot([1e-8,1e2],[1e-8,1e2],"--k")
axis([vx,Vx,vx,Vx])
xlabel("Predicted variance (theoretical data)")
ylabel("Measured variance")





