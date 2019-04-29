using DelimitedFiles,Statistics,SparseArrays,Dates

include("kuramoto.jl")

function sync(ntw::String,P0::Float64,M::Array{Float64,1},D::Array{Float64,1})
	@info "$(now()) -- Computing sync for P0 = $P0"
	
	Asp = readdlm(ntw*"_adj_mat.csv",',')
	
	n = round(Int,maximum(Asp[:,1]))
	m = Int(size(Asp)[1]/2)
	
	A = sparse(vec(Asp[:,1]),vec(Asp[:,2]),vec(Asp[:,3]))
	L = spdiagm(0 => vec(sum(A,dims=2))) - A
	B,w = L2B(L)
		
	P = P0*vec(readdlm("P_"*ntw*".csv",','))
	P .-= mean(P)
	
	x0 = zeros(2*n)
	
	x1,dx1 = kuramoto2_incidence(B,w,M,D,P,x0[1:n],x0[(n+1):(2*n)])
	
	writedlm("sync_states/"*ntw*"_sync_$P0.csv",x1,',')
	
	return x1,dx1
end



