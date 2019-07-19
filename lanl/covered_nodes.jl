using LightGraphs, SimpleWeightedGraphs

## Given a graph's sparse adjacency matrix (Asp) and a subset of its nodes (pmu_idx), return the list of nodes which appear in the shortest path between two of the selected nodes.

function covered_nodes(Asp::Array{Float64,2}, pmu_idx::Array{Int64,1})
	g = SimpleWeightedGraph(Array{Int64,1}(Asp[:,1]),Array{Int64,1}(Asp[:,2]),Asp[:,3])
	
	cn = []
	
	for i in 1:length(pmu_idx)-1
		for j in i+1:length(pmu_idx)
			cn = union(enumerate_paths(dijkstra_shortest_paths(g,pmu_idx[i]),pmu_idx[j]),cn)
		end
	end
	
	return cn
end







