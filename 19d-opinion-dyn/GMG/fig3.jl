using PyPlot, DelimitedFiles, Statistics

ϵs = vec(readdlm("add_data/fig3-erange.csv",','))

cmap = get_cmap("YlOrRd")

figure("fig3",(20,6))

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

subplot2grid((6,3),(0,0),rowspan=3)
PyPlot.plot(ϵs,vec(mean(p_pr,dims=1)),color=cmap(.9))
PyPlot.plot(ϵs,vec(mean(p_wta,dims=1)),color=cmap(.6))
PyPlot.plot(ϵs,vec(mean(p_sr,dims=1)),color=cmap(.3))
axis([0.,maximum(ϵs),0.025,0.145])
xticks(xticks()[1][1:end-1],[])
#xlabel("ϵ")
ylabel("effort (percentage of agents)")

div = p1[1,1]+p2[1,1]+p3[1,1]

subplot2grid((6,3),(3,0))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:]+p2[1,:];p3[1,:]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:]+p2[1,:]+p1[1,:];p3[1,:]+p2[1,:]]./div,color=cmap(.9))
ylabel("PR")
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])

subplot2grid((6,3),(4,0))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:]+p2[2,:];p3[2,:]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:]+p2[2,:]+p1[2,:];p3[2,:]+p2[2,:]]./div,color=cmap(.9))
ylabel("WTA")
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])

subplot2grid((6,3),(5,0))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p2[3,:];p3[3,:]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p2[3,:]+p1[3,:];p3[3,:]+p2[3,:]]./div,color=cmap(.9))
ylabel("SR")
axis([ϵs[1],ϵs[end],0.,1.])
xlabel("ϵ")

 #=
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
# =#

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
subplot2grid((6,3),(0,1),rowspan=3)
PyPlot.plot(ϵs,vec(mean(p_pr,dims=1)),color=cmap(.9))
PyPlot.plot(ϵs,vec(mean(p_wta,dims=1)),color=cmap(.6))
PyPlot.plot(ϵs,vec(mean(p_sr,dims=1)),color=cmap(.3))
axis([0.,maximum(ϵs),0.025,0.145])
xticks(xticks()[1][1:end-1],[])
#xlabel("ϵ")
ylabel("effort (percentage of agents)")

div = p1[1,1]+p2[1,1]+p3[1,1]

subplot2grid((6,3),(3,1))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:]+p2[1,:];p3[1,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:]+p2[1,:]+p1[1,:];p3[1,end:-1:1]+p2[1,end:-1:1]]./div,color=cmap(.9))
ylabel("PR")
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])

subplot2grid((6,3),(4,1))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:]+p2[2,:];p3[2,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:]+p2[2,:]+p1[2,:];p3[2,end:-1:1]+p2[2,end:-1:1]]./div,color=cmap(.9))
ylabel("WTA")
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])

subplot2grid((6,3),(5,1))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p2[3,:];p3[3,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p2[3,:]+p1[3,:];p3[3,end:-1:1]+p2[3,end:-1:1]]./div,color=cmap(.9))
ylabel("SR")
axis([ϵs[1],ϵs[end],0.,1.])
xlabel("ϵ")

 #=
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
# =#

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
subplot2grid((6,3),(0,2),rowspan=3)
PyPlot.plot(ϵs,vec(mean(p_pr,dims=1)),color=cmap(.9))
PyPlot.plot(ϵs,vec(mean(p_wta,dims=1)),color=cmap(.6))
PyPlot.plot(ϵs,vec(mean(p_sr,dims=1)),color=cmap(.3))
axis([0.,maximum(ϵs),0.025,0.145])
xticks(xticks()[1][1:end-1],[])
#xlabel("ϵ")
ylabel("effort (percentage of agents)")

div = p1[1,1]+p2[1,1]+p3[1,1]

subplot2grid((6,3),(3,2))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:]+p2[1,:];p3[1,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[1,:]+p2[1,:]+p1[1,:];p3[1,end:-1:1]+p2[1,end:-1:1]]./div,color=cmap(.9))
ylabel("PR")
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])

subplot2grid((6,3),(4,2))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:]+p2[2,:];p3[2,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[2,:]+p2[2,:]+p1[2,:];p3[2,end:-1:1]+p2[2,end:-1:1]]./div,color=cmap(.9))
ylabel("WTA")
axis([ϵs[1],ϵs[end],0.,1.])
xticks(xticks()[1][1:end-1],[])

subplot2grid((6,3),(5,2))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:];zeros(length(ϵs))]./div,color=cmap(.3))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p2[3,:];p3[3,end:-1:1]]./div,color=cmap(.6))
PyPlot.fill([ϵs;ϵs[end:-1:1]],[p3[3,:]+p2[3,:]+p1[3,:];p3[3,end:-1:1]+p2[3,end:-1:1]]./div,color=cmap(.9))
ylabel("SR")
axis([ϵs[1],ϵs[end],0.,1.])
xlabel("ϵ")

 #=
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
# =#

###################################################

l = 200
x = LinRange(-5,5,l)
figure("Distributions",(12,2))

p = .4
μ = 0.
Δ = 1.2
subplot(1,3,1)
y1 = p*exp.(-(x.+Δ.-μ).^2) + (1-p)*exp.(-(x.-Δ.-μ).^2)
PyPlot.fill([x;x[end:-1:1]],[y1;zeros(l)],color=[.7,.7,.7])
PyPlot.plot(x,y1,"--k")
PyPlot.plot([0,0],[-.1,1.],"k",linewidth=1)
PyPlot.plot([-5,5],[0,0],"k",linewidth=1)
axis([-5,5,-.1,.7])
xticks([])
yticks([])

p = .5
μ = 0.3
Δ = 1.2
subplot(1,3,2)
y1 = p*exp.(-(x.+Δ.-μ).^2) + (1-p)*exp.(-(x.-Δ.-μ).^2)
PyPlot.fill([x;x[end:-1:1]],[y1;zeros(l)],color=[.7,.7,.7])
PyPlot.plot(x,y1,"--k")
PyPlot.plot([0,0],[-.1,1.],"k",linewidth=1)
PyPlot.plot([-5,5],[0,0],"k",linewidth=1)
axis([-5,5,-.1,.7])
xticks([])
yticks([])

p = .5
μ = 0.3
Δ = 0.
subplot(1,3,3)
y1 = .6*(p*exp.(-(x.+Δ.-μ).^2) + (1-p)*exp.(-(x.-Δ.-μ).^2))
PyPlot.fill([x;x[end:-1:1]],[y1;zeros(l)],color=[.7,.7,.7])
PyPlot.plot(x,y1,"--k")
PyPlot.plot([0,0],[-.1,1.],"k",linewidth=1)
PyPlot.plot([-5,5],[0,0],"k",linewidth=1)
axis([-5,5,-.1,.7])
xticks([])
yticks([])
