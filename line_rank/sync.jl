using DelimitedFiles,Statistics,SparseArrays,Dates

include("kuramoto.jl")

function sync(ntw::String,P0::Float64,M::Array{Float64,1},D::Array{Float64,1},max_iter::Int64=50)
	@info "$(now()) -- Computing sync for P0 = $P0"
	
	Bsp = readdlm(ntw*"_data/"*ntw*"_inc_mat.csv",',')
	
	B = sparse(Bsp[:,1],Bsp[:,2],Bsp[:,3])
	n,m = size(B)
	
	w = Bsp[2*(1:m),4]
	
	P = P0*vec(readdlm(ntw*"_data/P_"*ntw*".csv",','))
	P .-= mean(P)
	
	x1,iter = NR_kuramoto(B,w,P)
	
	if iter >= max_iter
		@info "No sync state found..."
		writedlm("sync_states/"*ntw*"_sync_$P0.csv",["nope",],',')
	else
		writedlm("sync_states/"*ntw*"_sync_$P0.csv",x1,',')
	end
	
	return x1
end



