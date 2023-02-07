using PyPlot, Statistics, DelimitedFiles

cmap = get_cmap("YlOrRd")

figure("fig4sup",(20,6))

####################### p = 3 #################################################
subplot2grid((4,4),(0,0),colspan=1,rowspan=2)
ϵs = vec(readdlm("add_data/fig4a-erange.csv",','))
eff_extremist = vec(readdlm("add_data/fig4a-1-3-rm.csv",',') + readdlm("add_data/fig4a-3-3-rm.csv",','))
eff_centrist = vec(readdlm("add_data/fig4a-2-3-rm.csv",','))
eff_tot = eff_extremist+eff_centrist

PyPlot.plot(ϵs,eff_extremist,color=cmap(.9),"--",label="extremists")
PyPlot.plot(ϵs,eff_centrist,color=cmap(.3),"--",label="centrists")
PyPlot.plot(ϵs,eff_tot,"k")
axis([0.,maximum(ϵs),-.05*maximum(eff_tot),1.05*maximum(eff_tot)])
xticks(xticks()[1][1:end-1],[])
ylabel("ξ")
legend()
title("p = 3")

subplot2grid((4,4),(2,0),colspan=1,rowspan=1)
PP3 = 100*readdlm("add_data/fig4b-1st-3-rm.csv",',')'
sPP3 = zeros(1,length(ϵs))
for i in 1:3
	global sPP3 = [sPP3;sum(PP3[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP3[i,:];sPP3[i+1,151:-1:1]],color=cmap(i/3))
end
xticks(xticks()[1][1:end-1],[])
ylabel("percentage of wins")
axis([0.,maximum(ϵs),0,100])

subplot2grid((4,4),(3,0),colspan=1,rowspan=1)
PP3 = 100*readdlm("add_data/fig4b-2nd-3-rm.csv",',')'
sPP3 = zeros(1,length(ϵs))
for i in 1:3
	global sPP3 = [sPP3;sum(PP3[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP3[i,:];sPP3[i+1,151:-1:1]],color=cmap(i/3))
end
xlabel("ϵ")
ylabel("percentage of 2nd")
axis([0.,maximum(ϵs),0,100])


####################### p = 4 #################################################
subplot2grid((4,4),(0,1),colspan=1,rowspan=2)
eff_extremist = vec(readdlm("add_data/fig4a-1-4-rm.csv",',') + readdlm("add_data/fig4a-4-4-rm.csv",','))
eff_centrist = vec(readdlm("add_data/fig4a-2-4-rm.csv",',') + readdlm("add_data/fig4a-3-4-rm.csv",','))
eff_tot = eff_extremist+eff_centrist

PyPlot.plot(ϵs,eff_extremist,color=cmap(.9),"--",label="extremists")
PyPlot.plot(ϵs,eff_centrist,color=cmap(.3),"--",label="centrists")
PyPlot.plot(ϵs,eff_tot,"k")
axis([0.,maximum(ϵs),-.05*maximum(eff_tot),1.05*maximum(eff_tot)])
xticks(xticks()[1][1:end-1],[])
ylabel("ξ")
legend()
title("p = 4")

subplot2grid((4,4),(2,1),colspan=1,rowspan=1)
PP4 = 100*readdlm("add_data/fig4b-1st-4-rm.csv",',')'
sPP4 = zeros(1,length(ϵs))
for i in 1:4
	global sPP4 = [sPP4;sum(PP4[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP4[i,:];sPP4[i+1,151:-1:1]],color=cmap(i/4))
end
xticks(xticks()[1][1:end-1],[])
ylabel("percentage of wins")
axis([0.,maximum(ϵs),0,100])

subplot2grid((4,4),(3,1),colspan=1,rowspan=1)
PP4 = 100*readdlm("add_data/fig4b-2nd-4-rm.csv",',')'
sPP4 = zeros(1,length(ϵs))
for i in 1:4
	global sPP4 = [sPP4;sum(PP4[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP4[i,:];sPP4[i+1,151:-1:1]],color=cmap(i/4))
end
xlabel("ϵ")
ylabel("percentage of 2nd")
axis([0.,maximum(ϵs),0,100])


####################### p = 5 #################################################
subplot2grid((4,4),(0,2),colspan=1,rowspan=2)
eff_extremist = vec(readdlm("add_data/fig4a-1-5-rm.csv",',') + readdlm("add_data/fig4a-5-5-rm.csv",','))
eff_moderate = vec(readdlm("add_data/fig4a-2-5-rm.csv",',') + readdlm("add_data/fig4a-4-5-rm.csv",','))
eff_centrist = vec(readdlm("add_data/fig4a-3-5-rm.csv",','))
eff_tot = eff_extremist+eff_moderate+eff_centrist

PyPlot.plot(ϵs,eff_extremist,color=cmap(.9),"--",label="extremists")
PyPlot.plot(ϵs,eff_moderate,color=cmap(.6),"--",label="moderates")
PyPlot.plot(ϵs,eff_centrist,color=cmap(.3),"--",label="centrists")
PyPlot.plot(ϵs,eff_tot,"k")
axis([0.,maximum(ϵs),-.05*maximum(eff_tot),1.05*maximum(eff_tot)])
xticks(xticks()[1][1:end-1],[])
ylabel("ξ")
legend()
title("p = 5")
subplot2grid((4,4),(2,2),colspan=1,rowspan=1)
PP5 = 100*readdlm("add_data/fig4b-1st-5-rm.csv",',')'
sPP5 = zeros(1,length(ϵs))
for i in 1:5
	global sPP5 = [sPP5;sum(PP5[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP5[i,:];sPP5[i+1,151:-1:1]],color=cmap(i/5))
end
xticks(xticks()[1][1:end-1],[])
ylabel("percentage of wins")
axis([0.,maximum(ϵs),0,100])

subplot2grid((4,4),(3,2),colspan=1,rowspan=1)
PP5 = 100*readdlm("add_data/fig4b-2nd-5-rm.csv",',')'
sPP5 = zeros(1,length(ϵs))
for i in 1:5
	global sPP5 = [sPP5;sum(PP5[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP5[i,:];sPP5[i+1,151:-1:1]],color=cmap(i/5))
end
xlabel("ϵ")
ylabel("percentage of 2nd")
axis([0.,maximum(ϵs),0,100])


####################### p = 7 #################################################
subplot2grid((4,4),(0,3),colspan=1,rowspan=2)
eff_extremist = vec(readdlm("add_data/fig4a-1-7-rm.csv",',') + readdlm("add_data/fig4a-7-7-rm.csv",','))
eff_moderate = vec(readdlm("add_data/fig4a-2-7-rm.csv",',') + readdlm("add_data/fig4a-6-7-rm.csv",','))
eff_momoderate = vec(readdlm("add_data/fig4a-3-7-rm.csv",',') + readdlm("add_data/fig4a-5-7-rm.csv",','))
eff_centrist = vec(readdlm("add_data/fig4a-4-7-rm.csv",','))
eff_tot = eff_extremist+eff_moderate+eff_momoderate+eff_centrist

PyPlot.plot(ϵs,eff_extremist,color=cmap(.9),"--",label="extremists")
PyPlot.plot(ϵs,eff_moderate,color=cmap(.7),"--",label="moderates")
PyPlot.plot(ϵs,eff_momoderate,color=cmap(.5),"--",label="very moderates")
PyPlot.plot(ϵs,eff_centrist,color=cmap(.3),"--",label="centrists")
PyPlot.plot(ϵs,eff_tot,"k")
axis([0.,maximum(ϵs),-.05*maximum(eff_tot),1.05*maximum(eff_tot)])
xticks(xticks()[1][1:end-1],[])
ylabel("ξ")
legend()
title("p = 7")

subplot2grid((4,4),(2,3),colspan=1,rowspan=1)
PP7 = 100*readdlm("add_data/fig4b-1st-7-rm.csv",',')'
sPP7 = zeros(1,length(ϵs))
for i in 1:7
	global sPP7 = [sPP7;sum(PP7[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP7[i,:];sPP7[i+1,151:-1:1]],color=cmap(i/7))
end
xticks(xticks()[1][1:end-1],[])
ylabel("percentage of wins")
axis([0.,maximum(ϵs),0,100])

subplot2grid((4,4),(3,3),colspan=1,rowspan=1)
PP7 = 100*readdlm("add_data/fig4b-2nd-7-rm.csv",',')'
sPP7 = zeros(1,length(ϵs))
for i in 1:7
	global sPP7 = [sPP7;sum(PP7[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP7[i,:];sPP7[i+1,151:-1:1]],color=cmap(i/7))
end
xlabel("ϵ")
ylabel("percentage of 2nd")
axis([0.,maximum(ϵs),0,100])








figure("fig5sup",(20,6))

ϵs = vec(readdlm("add_data/fig4b-erange.csv",','))

####################### p = 3 #################################################
p = 3
p1 = readdlm("add_data/fig4b-n-1st-pr-wta-sr-$p-rm.csv",',')
p2 = readdlm("add_data/fig4b-n-2nd-pr-wta-sr-$p-rm.csv",',')
p3 = readdlm("add_data/fig4b-n-3rd-pr-wta-sr-$p-rm.csv",',')

div = p1[1,1]+p2[1,1]+p3[1,1]

subplot(3,4,1)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:]+p2[1,:];p3[1,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:]+p2[1,:]+p1[1,:];p3[1,end:-1:1]+p2[1,end:-1:1]]./div,color=cmap(.9))
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])
ylabel("PR")
title("p = 3")

subplot(3,4,5)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:]+p2[2,:];p3[2,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:]+p2[2,:]+p1[2,:];p3[2,end:-1:1]+p2[2,end:-1:1]]./div,color=cmap(.9))
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])
ylabel("WTA")

subplot(3,4,9)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p2[3,:];p3[3,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p2[3,:]+p1[3,:];p3[3,end:-1:1]+p2[3,end:-1:1]]./div,color=cmap(.9))
axis([ϵs[1],ϵs[end],0.,1.])
xlabel("ϵ")
ylabel("SR")


####################### p = 4 #################################################
p = 4
p1 = readdlm("add_data/fig4b-n-1st-pr-wta-sr-$p-rm.csv",',')
p2 = readdlm("add_data/fig4b-n-2nd-pr-wta-sr-$p-rm.csv",',')
p3 = readdlm("add_data/fig4b-n-3rd-pr-wta-sr-$p-rm.csv",',')

div = p1[1,1]+p2[1,1]+p3[1,1]

subplot(3,4,2)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:]+p2[1,:];p3[1,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:]+p2[1,:]+p1[1,:];p3[1,end:-1:1]+p2[1,end:-1:1]]./div,color=cmap(.9))
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])
ylabel("PR")
title("p = 4")

subplot(3,4,6)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:]+p2[2,:];p3[2,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:]+p2[2,:]+p1[2,:];p3[2,end:-1:1]+p2[2,end:-1:1]]./div,color=cmap(.9))
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])
ylabel("WTA")

subplot(3,4,10)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p2[3,:];p3[3,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p2[3,:]+p1[3,:];p3[3,end:-1:1]+p2[3,end:-1:1]]./div,color=cmap(.9))
axis([ϵs[1],ϵs[end],0.,1.])
xlabel("ϵ")
ylabel("SR")


####################### p = 5 #################################################
p = 5
p1 = readdlm("add_data/fig4b-n-1st-pr-wta-sr-$p-rm.csv",',')
p2 = readdlm("add_data/fig4b-n-2nd-pr-wta-sr-$p-rm.csv",',')
p3 = readdlm("add_data/fig4b-n-3rd-pr-wta-sr-$p-rm.csv",',')

div = p1[1,1]+p2[1,1]+p3[1,1]

subplot(3,4,3)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:]+p2[1,:];p3[1,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:]+p2[1,:]+p1[1,:];p3[1,end:-1:1]+p2[1,end:-1:1]]./div,color=cmap(.9))
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])
ylabel("PR")
title("p = 5")

subplot(3,4,7)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:]+p2[2,:];p3[2,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:]+p2[2,:]+p1[2,:];p3[2,end:-1:1]+p2[2,end:-1:1]]./div,color=cmap(.9))
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])
ylabel("WTA")

subplot(3,4,11)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p2[3,:];p3[3,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p2[3,:]+p1[3,:];p3[3,end:-1:1]+p2[3,end:-1:1]]./div,color=cmap(.9))
axis([ϵs[1],ϵs[end],0.,1.])
xlabel("ϵ")
ylabel("SR")


####################### p = 7 #################################################
p = 7
p1 = readdlm("add_data/fig4b-n-1st-pr-wta-sr-$p-rm.csv",',')
p2 = readdlm("add_data/fig4b-n-2nd-pr-wta-sr-$p-rm.csv",',')
p3 = readdlm("add_data/fig4b-n-3rd-pr-wta-sr-$p-rm.csv",',')

div = p1[1,1]+p2[1,1]+p3[1,1]

subplot(3,4,4)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:]+p2[1,:];p3[1,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:]+p2[1,:]+p1[1,:];p3[1,end:-1:1]+p2[1,end:-1:1]]./div,color=cmap(.9))
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])
ylabel("PR")
title("p = 7")

subplot(3,4,8)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:]+p2[2,:];p3[2,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:]+p2[2,:]+p1[2,:];p3[2,end:-1:1]+p2[2,end:-1:1]]./div,color=cmap(.9))
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])
ylabel("WTA")

subplot(3,4,12)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p2[3,:];p3[3,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p2[3,:]+p1[3,:];p3[3,end:-1:1]+p2[3,end:-1:1]]./div,color=cmap(.9))
axis([ϵs[1],ϵs[end],0.,1.])
xlabel("ϵ")
ylabel("SR")





