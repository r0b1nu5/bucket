using PyPlot, Statistics, DelimitedFiles

cmap = get_cmap("YlOrRd")

figure("fig4",(20,4.5))

subplot2grid((2,3),(0,0),colspan=1,rowspan=2)
ϵs = vec(readdlm("add_data/fig4a-erange.csv",','))
eff_extremist = vec(readdlm("add_data/fig4a-extremist.csv",','))
eff_moderate = vec(readdlm("add_data/fig4a-moderate.csv",','))
eff_centrist = vec(readdlm("add_data/fig4a-centrist.csv",','))
eff_tot = vec(readdlm("add_data/fig4a-total.csv",','))

PyPlot.plot(ϵs,eff_extremist,color=cmap(.9),"--",label="extremists")
PyPlot.plot(ϵs,eff_moderate,color=cmap(.6),"--",label="moderates")
PyPlot.plot(ϵs,eff_centrist,color=cmap(.3),"--",label="centrists")
PyPlot.plot(ϵs,eff_extremist+eff_moderate+eff_centrist,"k")
axis([0.,maximum(ϵs),-.1,2.5])
xlabel("ϵ")
ylabel("ξ")
legend()

subplot2grid((2,3),(0,1),colspan=1,rowspan=1)
PP6 = readdlm("add_data/fig4b-1st.csv",',')'
sPP6 = zeros(1,length(ϵs))
for i in 1:6
	global sPP6 = [sPP6;sum(PP6[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP6[i,:];sPP6[i+1,151:-1:1]],color=cmap(i/6))
end
xlabel("ϵ")
ylabel("percentage of wins")
axis([0.,maximum(ϵs),0,100])

subplot2grid((2,3),(1,1),colspan=1,rowspan=1)
PP6 = readdlm("add_data/fig4b-2nd.csv",',')'
sPP6 = zeros(1,length(ϵs))
for i in 1:6
	global sPP6 = [sPP6;sum(PP6[1:i,:],dims=1)]
	PyPlot.fill([ϵs;ϵs[end:-1:1]],[sPP6[i,:];sPP6[i+1,151:-1:1]],color=cmap(i/6))
end
xlabel("ϵ")
ylabel("percentage of 2nd")
axis([0.,maximum(ϵs),0,100])


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

