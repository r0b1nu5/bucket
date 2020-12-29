using DelimitedFiles, LightGraphs, SimpleWeightedGraphs, SparseArrays

# Construct a graph from the file containing the sparse Laplacian.
# "matrix": 	adj: uses the adjacency matrix
# 		lap: uses the Laplacian matrix
# "type":	simple: SimpleGraph
# 		directed: SimpleDiGraph
# 		weighted: SimpleWeightedGraph
# 		dweighted: SimpleWeightedDiGraph

function load_graph(file::String, matrix::String, type::String="simple")
	X = readdlm("file",',')
	I = Int.(X[:,1])
	J = Int.(X[:,2])
	V = X[:,3]
	l = length(I)

	if matrix == "adj"
		A = sparse(I,J,V)
	elseif matrix == "lap"
		ids = setdiff((1:l).*(I .!= J),[0.,])
		A = sparse(I[ids],J[ids],-V[ids])
	else
		@info "Wrong \"matrix\" argument."
		A = spzeros(1,1)
	end

	if type == "simple"
		return SimpleGraph(A)
	elseif type == "directed"
		return SimpleDiGraph(A)
	elseif type == "weighted"
		return SimpleWeightedGraph(A)
	elseif type == "dweighted"
		return SimpleWeightedDiGraph(A)
	else
		@info "Wrong \"type\" argument."
		return nothing
	end
end




