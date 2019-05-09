using PyPlot,DelimitedFiles

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
