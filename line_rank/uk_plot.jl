using PyPlot,DelimitedFiles,DelimitedFiles

include("../get_rgb.jl")

function uk_plot()
	@info("Plotting UK network...")
	
	coord = readdlm("uk_data/uk_grid_coord.csv",',')
	bord = readdlm("uk_data/uk_border_coord.csv",',')
	adj = Array{Int,2}(readdlm("uk_data/uk_adj_mat.csv",',')[:,1:2] .+ 1)
	
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
end


function uk_plot(dot_col::Array{Float64,1},line_col::Array{Float64,1}=[0.,0.,0.],bord_col::Array{Float64,1}=[0.,0.,0.])
	@info("Plotting UK network...")
	
	coord = readdlm("uk_data/uk_grid_coord.csv",',')
	bord = readdlm("uk_data/uk_border_coord.csv",',')
	adj = Array{Int,2}(readdlm("uk_data/uk_adj_mat.csv",',')[:,1:2] .+ 1)
	
	figure("UK",(8,10))
	PyPlot.plot(vec(bord[:,1]),vec(bord[:,2]),color=bord_col,linewidth=1)
	for i in 1:165
		id = 2*i-1
		PyPlot.plot([coord[adj[id,1],1],coord[adj[id,2],1]],[coord[adj[id,1],2],coord[adj[id,2],2]],color=line_col)
	end
	PyPlot.plot(vec(coord[:,1]),vec(coord[:,2]),"o",color=dot_col)
	for i in 1:120
		PyPlot.text(coord[i,1]+.005,coord[i,2]+.005,"$i")
	end
end

# iscolored[1]: "true" if dots are colored
# iscolored[2]: "true" if lines are colored
function uk_plot(iscolored::Array{Bool,1},dot_range::Array{Float64,1},line_range::Array{Float64,1})
	@info("Plotting UK network...")
	
	coord = readdlm("uk_data/uk_grid_coord.csv",',')
	bord = readdlm("uk_data/uk_border_coord.csv",',')
	adj = Array{Int,2}(readdlm("uk_data/uk_adj_mat.csv",',')[:,1:2] .+ 1)
	
	figure("UK",(8,10))
	PyPlot.plot(vec(bord[:,1]),vec(bord[:,2]),"k",linewidth=1)
	
	if iscolored[2]
		for i in 1:165
			id = 2*i-1
			PyPlot.plot([coord[adj[id,1],1],coord[adj[id,2],1]],[coord[adj[id,1],2],coord[adj[id,2],2]],color=get_rgb(line_range[i],minimum(line_range),maximum(line_range)))
		end
	else
		for i in 1:165
			id = 2*i-1
			PyPlot.plot([coord[adj[id,1],1],coord[adj[id,2],1]],[coord[adj[id,1],2],coord[adj[id,2],2]],"k")
		end
	end

	if iscolored[1]
		for i in 1:120
			PyPlot.plot(coord[i,1],coord[i,2],"o",color=get_rgb(dot_range[i],minimum(dot_range),maximum(dot_range)))
		end
	else
		PyPlot.plot(vec(coord[:,1]),vec(coord[:,2]),"ob")
	end

	for i in 1:120
		PyPlot.text(coord[i,1]+.005,coord[i,2]+.005,"$i")
	end
end

function uk_plot_ranked_lines(ranks::Array{Int64,1},ttl::String="UK, ?, ?, P0 = ?: Ranks of the lines",fig_size::Tuple{Int64,Int64}=(8,10))
	@info("Plotting UK network...")
	
	Bsp = readdlm("uk_data/uk_inc_mat.csv",',')
	B = sparse(Bsp[:,1],Bsp[:,2],Bsp[:,3])
	
	coord = readdlm("uk_data/uk_grid_coord.csv",',')
	bord = readdlm("uk_data/uk_border_coord.csv",',')
	
	rgb = get_rgb(1,length(ranks))
	
	lf = get_fignums()
	if length(lf) > 0
		num = maximum(lf)+1
	else
		num = 1
	end
	figure(num,fig_size)
	PyPlot.plot(vec(bord[:,1]),vec(bord[:,2]),"k",linewidth=1)
	for i in 1:length(ranks)
		idx = findnz(B[:,ranks[i]])[1]
		c1 = coord[idx[1],:]
		c2 = coord[idx[2],:]
		cm = (c1+c2)./2
		PyPlot.plot([c1[1],c2[1]],[c1[2],c2[2]],color=rgb[length(ranks)-i+1])
		PyPlot.plot(cm[1],cm[2],"s",markersize=10,color=rgb[length(ranks)-i+1])
		PyPlot.text(cm[1],cm[2],"$i",ha="center",va="center")
	end
	for i in setdiff(1:165,ranks)
		idx = findnz(B[:,i])[1]
		c1 = coord[idx[1],:]
		c2 = coord[idx[2],:]
		cm = (c1+c2)./2
		PyPlot.plot([c1[1],c2[1]],[c1[2],c2[2]],"k")
	end

	gen_idx = [14,33,36,44,51,58,95,104,112,119,120]
	con_idx = setdiff(1:120,gen_idx)	
	PyPlot.plot(vec(coord[con_idx,1]),vec(coord[con_idx,2]),"ok")
	PyPlot.plot(vec(coord[gen_idx,1]),vec(coord[gen_idx,2]),"sm")
	
	title(ttl)

end





