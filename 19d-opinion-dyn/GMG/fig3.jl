using PyPlot, DelimitedFiles, Statistics

ϵs = vec(readdlm("add_data/fig3-erange.csv",','))

cmap = get_cmap("YlOrRd")

figure("fig3")

###################################################

p_pr = readdlm("add_data/fig3a-pr.csv",',')
p_wta = readdlm("add_data/fig3a-wta.csv",',')
p_sr = readdlm("add_data/fig3a-sr.csv",',')

p_pr_1 = sum((p_pr .>= p_wta).*(p_pr .>= p_sr),dims=1)
p_wta_1 = sum((p_wta .> p_pr).*(p_wta .>= p_sr),dims=1)
p_sr_1 = sum((p_sr .> p_wta).*(p_sr .> p_pr),dims=1)
p1 = [p_pr_1;p_wta_1;p_sr_1]

p_pr_2 = sum((p_pr .>= p_wta).*(p_pr .< p_sr) + (p_pr .< p_wta).*(p_pr .>= p_sr),dims=1)
p_wta_2 = sum((p_wta .> p_pr).*(p_wta .< p_sr) + (p_wta .<= p_pr).*(p_wta .>= p_sr),dims=1)
p_sr_2 = sum((p_sr .> p_wta).*(p_sr .<= p_pr) + (p_sr .<= p_wta).*(p_sr .> p_pr),dims=1)
p2 = [p_pr_2;p_wta_2;p_sr_2]

p_pr_3 = sum((p_pr .< p_wta).*(p_pr .< p_sr),dims=1)
p_wta_3 = sum((p_wta .<= p_pr).*(p_wta .< p_sr),dims=1)
p_sr_3 = sum((p_sr .<= p_wta).*(p_sr .<= p_pr),dims=1)
p3 = [p_pr_3;p_wta_3;p_sr_3]

subplot(1,3,1)
PyPlot.plot(ϵs,vec(mean(p_pr,dims=1)),color=cmap(.9))
PyPlot.plot(ϵs,vec(mean(p_wta,dims=1)),color=cmap(.6))
PyPlot.plot(ϵs,vec(mean(p_sr,dims=1)),color=cmap(.3))
xlabel("ϵ")
ylabel("effort (percentage of agents)")

#matshow([p1;p2;p3],cmap=cmap,vmin=0.,vmax=1500.)
matshow([p1;p2;p3],cmap=cmap,vmin=0.,vmax=1100.)
colorbar()
PyPlot.plot([-.5,24.5],[2.5,2.5],"k")
PyPlot.plot([-.5,24.5],[5.5,5.5],"k")
yticks([0,1,2,3,4,5,6,7,8],["PR","WTA","SR","PR","WTA","SR","PR","WTA","SR"])
xlabel("ϵ")
###################################################

p_pr = readdlm("add_data/fig3b-pr.csv",',')
p_wta = readdlm("add_data/fig3b-wta.csv",',')
p_sr = readdlm("add_data/fig3b-sr.csv",',')

p_pr_1 = sum((p_pr .>= p_wta).*(p_pr .>= p_sr),dims=1)
p_wta_1 = sum((p_wta .> p_pr).*(p_wta .>= p_sr),dims=1)
p_sr_1 = sum((p_sr .> p_wta).*(p_sr .> p_pr),dims=1)
p1 = [p_pr_1;p_wta_1;p_sr_1]

p_pr_2 = sum((p_pr .>= p_wta).*(p_pr .< p_sr) + (p_pr .< p_wta).*(p_pr .>= p_sr),dims=1)
p_wta_2 = sum((p_wta .> p_pr).*(p_wta .< p_sr) + (p_wta .<= p_pr).*(p_wta .>= p_sr),dims=1)
p_sr_2 = sum((p_sr .> p_wta).*(p_sr .<= p_pr) + (p_sr .<= p_wta).*(p_sr .> p_pr),dims=1)
p2 = [p_pr_2;p_wta_2;p_sr_2]

p_pr_3 = sum((p_pr .< p_wta).*(p_pr .< p_sr),dims=1)
p_wta_3 = sum((p_wta .<= p_pr).*(p_wta .< p_sr),dims=1)
p_sr_3 = sum((p_sr .<= p_wta).*(p_sr .<= p_pr),dims=1)
p3 = [p_pr_3;p_wta_3;p_sr_3]

figure("fig3")
subplot(1,3,2)
PyPlot.plot(ϵs,vec(mean(p_pr,dims=1)),color=cmap(.9),label="prop. rep.")
PyPlot.plot(ϵs,vec(mean(p_wta,dims=1)),color=cmap(.6),label="w-t-a")
PyPlot.plot(ϵs,vec(mean(p_sr,dims=1)),color=cmap(.3),label="single rep.")
xlabel("ϵ")
ylabel("effort (percentage of agents)")
legend()

#matshow([p1;p2;p3],cmap=cmap,vmin=0.,vmax=1500.)
matshow([p1;p2;p3],cmap=cmap,vmin=0.,vmax=1100.)
colorbar()
PyPlot.plot([-.5,24.5],[2.5,2.5],"k")
PyPlot.plot([-.5,24.5],[5.5,5.5],"k")
yticks([0,1,2,3,4,5,6,7,8],["PR","WTA","SR","PR","WTA","SR","PR","WTA","SR"])
xlabel("ϵ")
###################################################

p_pr = readdlm("add_data/fig3c-pr.csv",',')
p_wta = readdlm("add_data/fig3c-wta.csv",',')
p_sr = readdlm("add_data/fig3c-sr.csv",',')

p_pr_1 = sum((p_pr .>= p_wta).*(p_pr .>= p_sr),dims=1)
p_wta_1 = sum((p_wta .> p_pr).*(p_wta .>= p_sr),dims=1)
p_sr_1 = sum((p_sr .> p_wta).*(p_sr .> p_pr),dims=1)
p1 = [p_pr_1;p_wta_1;p_sr_1]

p_pr_2 = sum((p_pr .>= p_wta).*(p_pr .< p_sr) + (p_pr .< p_wta).*(p_pr .>= p_sr),dims=1)
p_wta_2 = sum((p_wta .> p_pr).*(p_wta .< p_sr) + (p_wta .<= p_pr).*(p_wta .>= p_sr),dims=1)
p_sr_2 = sum((p_sr .> p_wta).*(p_sr .<= p_pr) + (p_sr .<= p_wta).*(p_sr .> p_pr),dims=1)
p2 = [p_pr_2;p_wta_2;p_sr_2]

p_pr_3 = sum((p_pr .< p_wta).*(p_pr .< p_sr),dims=1)
p_wta_3 = sum((p_wta .<= p_pr).*(p_wta .< p_sr),dims=1)
p_sr_3 = sum((p_sr .<= p_wta).*(p_sr .<= p_pr),dims=1)
p3 = [p_pr_3;p_wta_3;p_sr_3]

figure("fig3")
subplot(1,3,3)
PyPlot.plot(ϵs,vec(mean(p_pr,dims=1)),color=cmap(.9))
PyPlot.plot(ϵs,vec(mean(p_wta,dims=1)),color=cmap(.6))
PyPlot.plot(ϵs,vec(mean(p_sr,dims=1)),color=cmap(.3))
xlabel("ϵ")
ylabel("effort (percentage of agents)")

#matshow([p1;p2;p3],cmap=cmap,vmin=0.,vmax=1500.)
matshow([p1;p2;p3],cmap=cmap,vmin=0.,vmax=1100.)
colorbar()
PyPlot.plot([-.5,24.5],[2.5,2.5],"k")
PyPlot.plot([-.5,24.5],[5.5,5.5],"k")
yticks([0,1,2,3,4,5,6,7,8],["PR","WTA","SR","PR","WTA","SR","PR","WTA","SR"])
xlabel("ϵ")
###################################################


