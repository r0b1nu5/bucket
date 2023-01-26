using PyPlot, DelimitedFiles, Statistics

ϵs = vec(readdlm("add_data/fig3-erange.csv",','))

cmap = get_cmap("YlOrRd")

figure("fig3",(20,4.5))

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
axis([0.,maximum(ϵs),0.025,0.145])
xlabel("ϵ")
ylabel("effort (percentage of agents)")

#matshow([p1;p2;p3],cmap=cmap,vmin=0.,vmax=1500.)
matshow([p1;p2;p3],cmap=cmap,vmin=0.,vmax=1100.)
colorbar()
PyPlot.plot([-.5,24.5],[2.5,2.5],"k")
PyPlot.plot([-.5,24.5],[5.5,5.5],"k")
yticks([0,1,2,3,4,5,6,7,8],["PR","WTA","SR","PR","WTA","SR","PR","WTA","SR"])
xlabel("ϵ")

div = sum(p1[:,1])

figure("fig3b")
subplot(3,1,1)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p1[3,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p1[3,:]+p1[2,:];p1[3,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p1[3,:]+p1[2,:]+p1[1,:];p1[3,end:-1:1]+p1[2,end:-1:1]]./div,color=cmap(.9))
ylabel("Percentage of 1st position")
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])

subplot(3,1,2)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p2[3,:];zeros(length(ϵs))]./div,color=get_cmap("Blues")(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p2[3,:]+p2[2,:];p2[3,end:-1:1]]./div,color=get_cmap("Blues")(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p2[3,:]+p2[2,:]+p2[1,:];p2[3,end:-1:1]+p2[2,end:-1:1]]./div,color=get_cmap("Blues")(.9))
ylabel("Percentage of 2nd position")
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])

subplot(3,1,3)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:];zeros(length(ϵs))]./div,color=get_cmap("Greens")(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p3[2,:];p3[3,end:-1:1]]./div,color=get_cmap("Greens")(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p3[2,:]+p3[1,:];p3[3,end:-1:1]+p3[2,end:-1:1]]./div,color=get_cmap("Greens")(.9))
xlabel("ϵ")
ylabel("Percentage of 3rd position")
axis([ϵs[1],ϵs[end],0.,1.])

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
axis([0.,maximum(ϵs),0.025,0.145])
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

figure("fig3c")
subplot(3,1,1)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p1[3,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p1[3,:]+p1[2,:];p1[3,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p1[3,:]+p1[2,:]+p1[1,:];p1[3,end:-1:1]+p1[2,end:-1:1]]./div,color=cmap(.9))
ylabel("Percentage of 1st position")
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])

subplot(3,1,2)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p2[3,:];zeros(length(ϵs))]./div,color=get_cmap("Blues")(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p2[3,:]+p2[2,:];p2[3,end:-1:1]]./div,color=get_cmap("Blues")(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p2[3,:]+p2[2,:]+p2[1,:];p2[3,end:-1:1]+p2[2,end:-1:1]]./div,color=get_cmap("Blues")(.9))
ylabel("Percentage of 2nd position")
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])

subplot(3,1,3)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:];zeros(length(ϵs))]./div,color=get_cmap("Greens")(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p3[2,:];p3[3,end:-1:1]]./div,color=get_cmap("Greens")(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p3[2,:]+p3[1,:];p3[3,end:-1:1]+p3[2,end:-1:1]]./div,color=get_cmap("Greens")(.9))
xlabel("ϵ")
ylabel("Percentage of 3rd position")
axis([ϵs[1],ϵs[end],0.,1.])

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
axis([0.,maximum(ϵs),0.025,0.145])
xlabel("ϵ")
ylabel("effort (percentage of agents)")

#matshow([p1;p2;p3],cmap=cmap,vmin=0.,vmax=1500.)
matshow([p1;p2;p3],cmap=cmap,vmin=0.,vmax=1100.)
colorbar()
PyPlot.plot([-.5,24.5],[2.5,2.5],"k")
PyPlot.plot([-.5,24.5],[5.5,5.5],"k")
yticks([0,1,2,3,4,5,6,7,8],["PR","WTA","SR","PR","WTA","SR","PR","WTA","SR"])
xlabel("ϵ")

figure("fig3d")
subplot(3,1,1)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p1[3,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p1[3,:]+p1[2,:];p1[3,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p1[3,:]+p1[2,:]+p1[1,:];p1[3,end:-1:1]+p1[2,end:-1:1]]./div,color=cmap(.9))
ylabel("Percentage of 1st position")
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])

subplot(3,1,2)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p2[3,:];zeros(length(ϵs))]./div,color=get_cmap("Blues")(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p2[3,:]+p2[2,:];p2[3,end:-1:1]]./div,color=get_cmap("Blues")(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p2[3,:]+p2[2,:]+p2[1,:];p2[3,end:-1:1]+p2[2,end:-1:1]]./div,color=get_cmap("Blues")(.9))
ylabel("Percentage of 2nd position")
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])

subplot(3,1,3)
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:];zeros(length(ϵs))]./div,color=get_cmap("Greens")(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p3[2,:];p3[3,end:-1:1]]./div,color=get_cmap("Greens")(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p3[2,:]+p3[1,:];p3[3,end:-1:1]+p3[2,end:-1:1]]./div,color=get_cmap("Greens")(.9))
xlabel("ϵ")
ylabel("Percentage of 3rd position")
axis([ϵs[1],ϵs[end],0.,1.])


###################################################



