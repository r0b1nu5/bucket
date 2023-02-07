using PyPlot, Statistics, DelimitedFiles

cmap = get_cmap("YlOrRd")

ϵs = vec(readdlm("add_data/fig4a-erange.csv",','))

figure("fig4sup",(20,15))

# 3 parties
subplot2grid((6,4),(0,0),colspan=1,rowspan=2)
eff_extremist = vec(readdlm("add_data/fig4sup-extremist-3.csv",','))
eff_centrist = vec(readdlm("add_data/fig4sup-centrist-3.csv",','))

PyPlot.plot(ϵs,eff_extremist+eff_centrist,"k")
PyPlot.plot(ϵs,eff_extremist,color=cmap(.9),"--",label="extremists")
PyPlot.plot(ϵs,eff_centrist,color=cmap(.3),"--",label="centrists")
axis([0.,maximum(ϵs),-.1,22])
xlabel("ϵ")
ylabel("ξ")
legend()

# 4 parties
subplot2grid((6,4),(0,1),colspan=1,rowspan=2)
eff_extremist = vec(readdlm("add_data/fig4sup-extremist-4.csv",','))
eff_centrist = vec(readdlm("add_data/fig4sup-centrist-4.csv",','))

PyPlot.plot(ϵs,eff_extremist+eff_centrist,"k")
PyPlot.plot(ϵs,eff_extremist,color=cmap(.9),"--",label="extremists")
PyPlot.plot(ϵs,eff_centrist,color=cmap(.3),"--",label="centrists")
axis([0.,maximum(ϵs),-.1,2.5])
xlabel("ϵ")
ylabel("ξ")
legend()

# 5 parties
subplot2grid((6,4),(0,2),colspan=1,rowspan=2)
eff_extremist = vec(readdlm("add_data/fig4sup-extremist-5.csv",','))
eff_moderate = vec(readdlm("add_data/fig4sup-moderate-5.csv",','))
eff_centrist = vec(readdlm("add_data/fig4sup-centrist-5.csv",','))

PyPlot.plot(ϵs,eff_extremist+eff_centrist+eff_moderate,"k")
PyPlot.plot(ϵs,eff_extremist,color=cmap(.9),"--",label="extremists")
PyPlot.plot(ϵs,eff_moderate,color=cmap(.6),"--",label="moderates")
PyPlot.plot(ϵs,eff_centrist,color=cmap(.3),"--",label="centrists")
axis([0.,maximum(ϵs),-.1,2.8])
xlabel("ϵ")
ylabel("ξ")
legend()

# 7 parties
subplot2grid((6,4),(0,3),colspan=1,rowspan=2)
eff_extremist = vec(readdlm("add_data/fig4sup-extremist-7.csv",','))
eff_moderate = vec(readdlm("add_data/fig4sup-moderate-7.csv",','))
eff_momoderate = vec(readdlm("add_data/fig4sup-momoderate-7.csv",','))
eff_centrist = vec(readdlm("add_data/fig4sup-centrist-7.csv",','))

PyPlot.plot(ϵs,eff_extremist+eff_centrist+eff_moderate+eff_momoderate,"k")
PyPlot.plot(ϵs,eff_extremist,color=cmap(.9),"--",label="extremists")
PyPlot.plot(ϵs,eff_moderate,color=cmap(.7),"--",label="moderates")
PyPlot.plot(ϵs,eff_momoderate,color=cmap(.5),"--",label="very moderates")
PyPlot.plot(ϵs,eff_centrist,color=cmap(.3),"--",label="centrists")
axis([0.,maximum(ϵs),-.1,3])
xlabel("ϵ")
ylabel("ξ")
legend()

##########################################################################################

# 3 parties
subplot2grid((6,4),(2,0),colspan=1,rowspan=1)
PP3 = readdlm("add_data/fig4sup-1st-3.csv",',')'
sPP3 = zeros(1,length(ϵs))
for i in 1:3
	global sPP3 = [sPP3;sum(PP3[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP3[i,:];sPP3[i+1,151:-1:1]],color=cmap(i/3))
end
xlabel("ϵ")
ylabel("percentage of wins")
axis([0.,maximum(ϵs),0,100])

subplot2grid((6,4),(3,0),colspan=1,rowspan=1)
PP3 = readdlm("add_data/fig4sup-2nd-3.csv",',')'
sPP3 = zeros(1,length(ϵs))
for i in 1:3
	global sPP3 = [sPP3;sum(PP3[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP3[i,:];sPP3[i+1,151:-1:1]],color=cmap(i/3))
end
xlabel("ϵ")
ylabel("percentage of 2nd")
axis([0.,maximum(ϵs),0,100])

# 4 parties
subplot2grid((6,4),(2,1),colspan=1,rowspan=1)
PP4 = readdlm("add_data/fig4sup-1st-4.csv",',')'
sPP4 = zeros(1,length(ϵs))
for i in 1:4
	global sPP4 = [sPP4;sum(PP4[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP4[i,:];sPP4[i+1,151:-1:1]],color=cmap(i/4))
end
xlabel("ϵ")
ylabel("percentage of wins")
axis([0.,maximum(ϵs),0,100])

subplot2grid((6,4),(3,1),colspan=1,rowspan=1)
PP4 = readdlm("add_data/fig4sup-2nd-4.csv",',')'
sPP4 = zeros(1,length(ϵs))
for i in 1:4
	global sPP4 = [sPP4;sum(PP4[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP4[i,:];sPP4[i+1,151:-1:1]],color=cmap(i/4))
end
xlabel("ϵ")
ylabel("percentage of 2nd")
axis([0.,maximum(ϵs),0,100])

# 5 parties
subplot2grid((6,4),(2,2),colspan=1,rowspan=1)
PP5 = readdlm("add_data/fig4sup-1st-5.csv",',')'
sPP5 = zeros(1,length(ϵs))
for i in 1:5
	global sPP5 = [sPP5;sum(PP5[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP5[i,:];sPP5[i+1,151:-1:1]],color=cmap(i/5))
end
xlabel("ϵ")
ylabel("percentage of wins")
axis([0.,maximum(ϵs),0,100])

subplot2grid((6,4),(3,2),colspan=1,rowspan=1)
PP5 = readdlm("add_data/fig4sup-2nd-5.csv",',')'
sPP5 = zeros(1,length(ϵs))
for i in 1:5
	global sPP5 = [sPP5;sum(PP5[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP5[i,:];sPP5[i+1,151:-1:1]],color=cmap(i/5))
end
xlabel("ϵ")
ylabel("percentage of 2nd")
axis([0.,maximum(ϵs),0,100])

# 7 parties
subplot2grid((6,4),(2,3),colspan=1,rowspan=1)
PP7 = readdlm("add_data/fig4sup-1st-7.csv",',')'
sPP7 = zeros(1,length(ϵs))
for i in 1:7
	global sPP7 = [sPP7;sum(PP7[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP7[i,:];sPP7[i+1,151:-1:1]],color=cmap(i/7))
end
xlabel("ϵ")
ylabel("percentage of wins")
axis([0.,maximum(ϵs),0,100])

subplot2grid((6,4),(3,3),colspan=1,rowspan=1)
PP7 = readdlm("add_data/fig4sup-2nd-7.csv",',')'
sPP7 = zeros(1,length(ϵs))
for i in 1:7
	global sPP7 = [sPP7;sum(PP7[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP7[i,:];sPP7[i+1,151:-1:1]],color=cmap(i/7))
end
xlabel("ϵ")
ylabel("percentage of 2nd")
axis([0.,maximum(ϵs),0,100])



# 3 parties
eff_sr = vec(readdlm("add_data/fig4sup-sr-3.csv",','))
eff_wta = vec(readdlm("add_data/fig4sup-wta-3.csv",','))
eff_pr = vec(readdlm("add_data/fig4sup-pr-3.csv",','))

subplot2grid((6,4),(4,0),colspan=1,rowspan=2)
PyPlot.plot(ϵs,eff_sr,color=cmap(.3),label="single rep.")
PyPlot.plot(ϵs,eff_wta,color=cmap(.6),label="winner takes all")
PyPlot.plot(ϵs,eff_pr,color=cmap(.9),label="prop. rep.")
axis([0.,maximum(ϵs),.1,4])
xlabel("ϵ")
ylabel("ξ")
legend()

# 4 parties
eff_sr = vec(readdlm("add_data/fig4sup-sr-4.csv",','))
eff_wta = vec(readdlm("add_data/fig4sup-wta-4.csv",','))
eff_pr = vec(readdlm("add_data/fig4sup-pr-4.csv",','))

subplot2grid((6,4),(4,1),colspan=1,rowspan=2)
PyPlot.plot(ϵs,eff_sr,color=cmap(.3),label="single rep.")
PyPlot.plot(ϵs,eff_wta,color=cmap(.6),label="winner takes all")
PyPlot.plot(ϵs,eff_pr,color=cmap(.9),label="prop. rep.")
axis([0.,maximum(ϵs),.2,1.5])
xlabel("ϵ")
ylabel("ξ")
legend()

# 5 parties
eff_sr = vec(readdlm("add_data/fig4sup-sr-5.csv",','))
eff_wta = vec(readdlm("add_data/fig4sup-wta-5.csv",','))
eff_pr = vec(readdlm("add_data/fig4sup-pr-5.csv",','))

subplot2grid((6,4),(4,2),colspan=1,rowspan=2)
PyPlot.plot(ϵs,eff_sr,color=cmap(.3),label="single rep.")
PyPlot.plot(ϵs,eff_wta,color=cmap(.6),label="winner takes all")
PyPlot.plot(ϵs,eff_pr,color=cmap(.9),label="prop. rep.")
axis([0.,maximum(ϵs),.2,1.6])
xlabel("ϵ")
ylabel("ξ")
legend()

# 7 parties
eff_sr = vec(readdlm("add_data/fig4sup-sr-7.csv",','))
eff_wta = vec(readdlm("add_data/fig4sup-wta-7.csv",','))
eff_pr = vec(readdlm("add_data/fig4sup-pr-7.csv",','))

subplot2grid((6,4),(4,3),colspan=1,rowspan=2)
PyPlot.plot(ϵs,eff_sr,color=cmap(.3),label="single rep.")
PyPlot.plot(ϵs,eff_wta,color=cmap(.6),label="winner takes all")
PyPlot.plot(ϵs,eff_pr,color=cmap(.9),label="prop. rep.")
axis([0.,maximum(ϵs),.2,1.7])
xlabel("ϵ")
ylabel("ξ")
legend()


