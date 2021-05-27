using PyPlot, DelimitedFiles, SparseArrays

# No weights
function show_graph(ntw::String, base_file::String="./ntw_data/")
	x,y,L = load_graph_data(ntw,base_file)

	n = length(x)
	I,J,V = findnz(L)
	ids = setdiff((I .< J).*(1:length(I)),[0,])

	figure(ntw)
	for k in 1:length(ids)
		i = I[ids[k]]
		j = J[ids[k]]
		
		PyPlot.plot([x[i],x[j]],[y[i],y[j]],"-k")
	end
	for i in 1:n
		PyPlot.plot(x[i],y[i],"ok")
	end

	return nothing
end

# Weights on the nodes
# Weights are assumed to be rates, i.e., w>1 means "overloaded"
function show_graph(ntw::String, nw::Array{Float64,1}, base_file::String="./ntw_data/")
	x,y,L = load_graph_data(ntw,base_file)

	n = length(x)
	I,J,V = findnz(L)
	ids = setdiff((I .< J).*(1:length(I)),[0,])

	cn = (nw .> .5) + (nw .> .8) + (nw .> 1.) .+ 1
	col = ["C2","C0","C1","C3"]

	figure(ntw)
	for k in 1:length(ids)
		i = I[ids[k]]
		j = J[ids[k]]
		
		PyPlot.plot([x[i],x[j]],[y[i],y[j]],"-k")
	end
	for i in 1:n
		PyPlot.plot(x[i],y[i],"o",color=col[cn[i]])
	end

	return nothing
end

# Weights on the edges
# Weights are assumed to be rates, i.e., w>1 means "overloaded"
function show_graph(ntw::String, lw::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int64}}, base_file::String="./ntw_data/")
	x,y,L = load_graph_data(ntw,base_file)

	n = length(x)
	I,J,V = findnz(L)
	ids = setdiff((I .< J).*(1:length(I)),[0,])

	cl = (lw .> .5) + (lw .> .8) + (lw .> 1.) .+ 1
	col = ["C2","C0","C1","C3"]

	figure(ntw)
	for k in 1:length(ids)
		i = I[ids[k]]
		j = J[ids[k]]
		
		PyPlot.plot([x[i],x[j]],[y[i],y[j]],color=col[cl[i,j]])
	end
	for i in 1:n
		PyPlot.plot(x[i],y[i],"ok")
	end

	return nothing
end

# Weights on the nodes and edges
# Weights are assumed to be rates, i.e., w>1 means "overloaded"
function show_graph(ntw::String, nw::Array{Float64,1}, lw::Union{Array{Float64,2},SparseMatrixCSC{Float64,Int64}}, base_file::String="./ntw_data/")
	x,y,L = load_graph_data(ntw,base_file)

	n = length(x)
	I,J,V = findnz(L)
	ids = setdiff((I .< J).*(1:length(I)),[0,])

	cn = (nw .> .5) + (nw .> .8) + (nw .> 1.) .+ 1
	cl = (lw .> .5) + (lw .> .8) + (lw .> 1.) .+ 1
	col = ["C2","C0","C1","C3"]

	figure(ntw)
	for k in 1:length(ids)
		i = I[ids[k]]
		j = J[ids[k]]
		
		PyPlot.plot([x[i],x[j]],[y[i],y[j]],color=col[cl[i,j]])
	end
	for i in 1:n
		PyPlot.plot(x[i],y[i],"o",color=col[cn[i]])
	end

	return nothing
end




function load_graph_data(ntw::String, base_file::String="./ntw_data/")
	xy = readdlm(base_file*ntw*"_xy.csv",',')
	x = xy[:,1]
	y = xy[:,2]

	Lsp = readdlm(base_file*ntw*"_Lsp.csv",',')
	I = Int64.(Lsp[:,1])
	J = Int64.(Lsp[:,2])
	V = Lsp[:,3]

	return x,y,sparse(I,J,V)
end


