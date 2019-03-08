using PyPlot,DelimitedFiles

@info("Plotting UK network...")

coord = readdlm("uk_grid_coord.csv",',')
bord = readdlm("uk_border_coord.csv",',')
adj = Array{Int,2}(readdlm("uk_adj_mat.csv",',')[:,1:2] .+ 1)

figure("UK",(8,10))
PyPlot.plot(vec(bord[:,1]),vec(bord[:,2]),"k",linewidth=1)
for i in 1:165
	id = 2*i-1
	PyPlot.plot([coord[adj[id,1],1],coord[adj[id,2],1]],[coord[adj[id,1],2],coord[adj[id,2],2]],"k")
end
PyPlot.plot(vec(coord[:,1]),vec(coord[:,2]),"ob")
for i in 1:120
	PyPlot.text(coord[i,1]+.005,coord[i,2]+.005,"$i")
end

