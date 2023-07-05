using PyPlot, Statistics, DelimitedFiles

cmap = get_cmap("YlOrRd")

figure("fig4",(6,9))

subplot2grid((4,1),(0,0),colspan=1,rowspan=2)
ϵs = vec(readdlm("add_data/fig4a-erange.csv",','))[1:121]
#eff_extremist = vec(readdlm("add_data/fig4a-extremist.csv",','))
eff_extremist = 100*vec(readdlm("add_data/fig4a-1-6-rm.csv",',') + readdlm("add_data/fig4a-6-6-rm.csv",','))[1:121]
#eff_moderate = vec(readdlm("add_data/fig4a-moderate.csv",','))
eff_moderate = 100*vec(readdlm("add_data/fig4a-2-6-rm.csv",',') + readdlm("add_data/fig4a-5-6-rm.csv",','))[1:121]
#eff_centrist = vec(readdlm("add_data/fig4a-centrist.csv",','))
eff_centrist = 100*vec(readdlm("add_data/fig4a-3-6-rm.csv",',') + readdlm("add_data/fig4a-4-6-rm.csv",','))[1:121]
eff_tot = vec(readdlm("add_data/fig4a-total.csv",','))

PyPlot.plot(ϵs,eff_extremist,color=cmap(.9),"--",label="extremists")
PyPlot.plot(ϵs,eff_moderate,color=cmap(.6),"--",label="moderates")
PyPlot.plot(ϵs,eff_centrist,color=cmap(.3),"--",label="centrists")
PyPlot.plot(ϵs,eff_extremist+eff_moderate+eff_centrist,"k",label="total")
axis([0.,maximum(ϵs),-.5,12])
xticks(xticks()[1][1:end-1],[])
ylabel("ξ")
legend()

subplot2grid((4,1),(2,0),colspan=1,rowspan=1)
#PP6 = readdlm("add_data/fig4b-1st.csv",',')'
PP6 = 100*readdlm("add_data/fig4b-1st-6-rm.csv",',')'[:,1:121]
sPP6 = zeros(1,length(ϵs))
for i in 1:6
	global sPP6 = [sPP6;sum(PP6[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[121:-1:1]],[sPP6[i,:];sPP6[i+1,121:-1:1]],color=cmap(i/6))
end
xticks(xticks()[1][1:end-1],[])
ylabel("percentage of wins")
axis([0.,maximum(ϵs),0,100])

subplot2grid((4,1),(3,0),colspan=1,rowspan=1)
#PP6 = readdlm("add_data/fig4b-2nd.csv",',')'
PP6 = 100*readdlm("add_data/fig4b-2nd-6-rm.csv",',')'[:,1:121]
sPP6 = zeros(1,length(ϵs))
for i in 1:6
	global sPP6 = [sPP6;sum(PP6[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP6[i,:];sPP6[i+1,121:-1:1]],color=cmap(i/6))
end
xlabel("ϵ")
ylabel("percentage of 2nd")
axis([0.,maximum(ϵs),0,100])


# #=
nep = 30
n1n2 = readdlm("add_data/n1n2-$p-$nep.csv",',')
n1n2n = n1n2./repeat(sum(n1n2,dims=1),p,1)
matshow(n1n2,cmap=cmap)
colorbar()
matshow(n1n2n,cmap=cmap)
colorbar()
# =#
 #=
subplot2grid((2,3),(0,2),colspan=1,rowspan=2)
eff_sr = vec(readdlm("add_data/fig4c-sr.csv",','))
eff_wta = vec(readdlm("add_data/fig4c-wta.csv",','))
eff_pr = vec(readdlm("add_data/fig4c-pr.csv",','))
PyPlot.plot(ϵs,eff_sr,color=cmap(.3),label="single rep.")
PyPlot.plot(ϵs,eff_wta,color=cmap(.6),label="winner takes all")
PyPlot.plot(ϵs,eff_pr,color=cmap(.9),label="prop. rep.")
axis([0.,maximum(ϵs),.2,1.8])
xlabel("ϵ")
ylabel("ξ")
legend()
# =#

figure("fig5")
ϵs = vec(readdlm("add_data/fig4b-erange.csv",','))[1:41]
p1 = readdlm("add_data/fig4b-n-1st-pr-wta-sr-6-rm.csv",',')[:,1:41]
p2 = readdlm("add_data/fig4b-n-2nd-pr-wta-sr-6-rm.csv",',')[:,1:41]
p3 = readdlm("add_data/fig4b-n-3rd-pr-wta-sr-6-rm.csv",',')[:,1:41]

div = p1[1,1]+p2[1,1]+p3[1,1]

subplot(3,1,1)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:]+p2[1,:];p3[1,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:]+p2[1,:]+p1[1,:];p3[1,end:-1:1]+p2[1,end:-1:1]]./div,color=cmap(.9))
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])
yticks([0.,.25,.5,.75,1.],["0.0","","0.5","","1.0"])
ylabel("PR")

subplot(3,1,2)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:]+p2[2,:];p3[2,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:]+p2[2,:]+p1[2,:];p3[2,end:-1:1]+p2[2,end:-1:1]]./div,color=cmap(.9))
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])
yticks([0.,.25,.5,.75,1.],["0.0","","0.5","","1.0"])
ylabel("WTA")

subplot(3,1,3)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p2[3,:];p3[3,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p2[3,:]+p1[3,:];p3[3,end:-1:1]+p2[3,end:-1:1]]./div,color=cmap(.9))
axis([ϵs[1],ϵs[end],0.,1.])
yticks([0.,.25,.5,.75,1.],["0.0","","0.5","","1.0"])
xlabel("ϵ")
ylabel("SR")


