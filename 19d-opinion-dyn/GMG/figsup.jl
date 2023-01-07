using DelimitedFiles, PyPlot, Statistics

influence = vec(readdlm("add_data/figsup-influence.csv",','))[1:17]

r_mean_d1 = vec(readdlm("add_data/figsup-corrmatmean-D1.csv",','))[1:17]
r_std_d1 = vec(readdlm("add_data/figsup-corrmatstd-D1.csv",','))[1:17]
r_max_d1 = vec(readdlm("add_data/figsup-corrmatmax-D1.csv",','))[1:17]

r_mean_d2 = vec(readdlm("add_data/figsup-corrmatmean-D2.csv",','))[1:17]
r_std_d2 = vec(readdlm("add_data/figsup-corrmatstd-D2.csv",','))[1:17]
r_max_d2 = vec(readdlm("add_data/figsup-corrmatmax-D2.csv",','))[1:17]

r_mean_d3 = vec(readdlm("add_data/figsup-corrmatmean-D3.csv",','))[1:17]
r_std_d3 = vec(readdlm("add_data/figsup-corrmatstd-D3.csv",','))[1:17]
r_max_d3 = vec(readdlm("add_data/figsup-corrmatmax-D3.csv",','))[1:17]

cmap = get_cmap("YlOrRd")

figure("figsup",(8,4))

PyPlot.plot(influence,r_mean_d1,"o",color=cmap(.3),label="D1")
for i in 1:length(influence)
	x = influence[i]
	y = r_mean_d1[i]
	z = r_std_d1[i]
	PyPlot.plot([x,x],[y-z,y+z],color=cmap(.3))
end
m,id = findmax(r_mean_d1)
PyPlot.plot([influence[id],influence[id]],[.35,.67],"--",color=cmap(.3))

PyPlot.plot(influence.+.1,r_mean_d2,"o",color=cmap(.6),label="D2")
for i in 1:length(influence)
	x = influence[i]
	y = r_mean_d2[i]
	z = r_std_d2[i]
	PyPlot.plot([x,x].+.1,[y-z,y+z],color=cmap(.6))
end
m,id = findmax(r_mean_d2)
PyPlot.plot([influence[id],influence[id]].+.1,[.35,.67],"--",color=cmap(.6))

PyPlot.plot(influence.+.2,r_mean_d3,"o",color=cmap(.9),label="D3")
for i in 1:length(influence)
	x = influence[i]
	y = r_mean_d3[i]
	z = r_std_d3[i]
	PyPlot.plot([x,x].+.2,[y-z,y+z],color=cmap(.9))
end
m,id = findmax(r_mean_d3)
PyPlot.plot([influence[id],influence[id]].+.2,[.35,.67],"--",color=cmap(.9))

axis([-1,26,.35,.67])
xlabel("Influence budget [%]")
ylabel("Pearson correlation")
legend()






