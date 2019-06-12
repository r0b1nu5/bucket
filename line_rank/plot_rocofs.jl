using PyPlot,DelimitedFiles,Statistics


 #=
Es = zeros(3,170)
vars = zeros(3,170)
EE = Array{Array{Float64,2},1}()
vv = Array{Array{Float64,2},1}()
for i in 1:9
	global EE,vv
	push!(EE,readdlm("data/Es_$i.csv",','))
	push!(vv,readdlm("data/Es_$i.csv",','))
end
for i in 1:9
	global Es,EE
	Es += EE[i]./9
end
vars[1:2,:] = vv[1][1:2,:]
for k in 1:170
	global vars,vv,EE,Es
	for i in 1:9
		vars[3,k] += 100/900*(vv[i][3,k] + (EE[i][3,k]-Es[3,k])^2)
	end
end
# =#

for n in [2206] # [102,103,104,105]	
 ##=
Es = readdlm("data/Es_$n.csv",',')
vars = readdlm("data/vars_$n.csv",',')
# =#

figure()
subplot(131)
PyPlot.plot(abs.(Es[1,:]),abs.(Es[3,:]),"x",label="Node i (r = $(round(cor(Es[1,:],Es[3,:]),digits=4)))",color=def_col[1])
PyPlot.plot([1e-8,100],[1e-8,100],"--k")
xlabel("Predicted expectation")
ylabel("Measured expectation")
legend()
axis([minimum(abs.(Es)),maximum(abs.(Es)),minimum(abs.(Es)),maximum(abs.(Es))])
subplot(132)
PyPlot.loglog(abs.(Es[1,:]),abs.(Es[3,:]),"x",label="Node i (r = $(round(cor(Es[1,:],Es[3,:]),digits=4)))",color=def_col[1])
PyPlot.plot([1e-8,100],[1e-8,100],"--k")
xlabel("Predicted expectation")
ylabel("Measured expectation")
legend()
axis([minimum(abs.(Es)),maximum(abs.(Es)),minimum(abs.(Es)),maximum(abs.(Es))])
subplot(133)
#=
PyPlot.plot(abs.(Es[2,:]),abs.(Es[3,:]),"x",label="Node j (r = $(round(cor(Es[2,:],Es[3,:]),digits=4)))",color=def_col[2])
PyPlot.plot([1e-8,100],[1e-8,100],"--k")
xlabel("Predicted expectation")
ylabel("Measured expectation")
legend()
axis([minimum(abs.(Es)),maximum(abs.(Es)),minimum(abs.(Es)),maximum(abs.(Es))])

figure()
=#
PyPlot.semilogy(sort(abs.(Es[1,:]-Es[3,:])))
ylabel("Error")
xlabel("Line rank")

figure()
subplot(131)
PyPlot.plot(vars[1,:],vars[3,:],"x",label="Node i (r = $(round(cor(vars[1,:],vars[3,:]),digits=4)))",color=def_col[2])
PyPlot.plot((10.).^(-8:2),(10.).^(-8:2),"--k")
xlabel("Predicted variance")
ylabel("Measured variance")
legend()
axis([minimum(vars),maximum(vars),minimum(vars),maximum(vars)])
subplot(132)
PyPlot.loglog(vars[1,:],vars[3,:],"x",label="Node i (r = $(round(cor(vars[1,:],vars[3,:]),digits=4)))",color=def_col[2])
PyPlot.plot((10.).^(-8:2),(10.).^(-8:2),"--k")
xlabel("Predicted variance")
ylabel("Measured variance")
legend()
axis([minimum(vars),maximum(vars),minimum(vars),maximum(vars)])
subplot(133)
#=
PyPlot.plot(vars[2,:],vars[3,:],"x",label="Node j (r = $(round(cor(vars[2,:],vars[3,:]),digits=4)))",color=def_col[2])
PyPlot.plot((10.).^(-8:2),(10.).^(-8:2),"--k")
xlabel("Predicted variance")
ylabel("Measured variance")
legend()
axis([minimum(vars),maximum(vars),minimum(vars),maximum(vars)])

figure()
=#
PyPlot.semilogy(sort(abs.(vars[1,:]-vars[3,:])),color=def_col[2])
ylabel("Error")
xlabel("Line rank")
 #= 
figure("all 1")
subplot(131)
PyPlot.plot(Es[1,:],"o",color=def_col[n-101])
subplot(132)
PyPlot.plot(Es[2,:],"o",color=def_col[n-101])
subplot(133)
PyPlot.plot(Es[3,:],"o",color=def_col[n-101])
# =#

end






