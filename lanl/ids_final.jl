using DelimitedFiles

include("load_pen.jl")

for date in dates
	Xs = load_pen(date)
	n = Int(size(Xs)[1]/2)
	
	Xs = Xs[[1:n-1;n+1:2*n-1],:]
	n = n-1
	
	ii = Array{Int64,1}()
	for i in 1:n-1
		for j in i+1:n
			if Xs[i,:] == Xs[j,:] || Xs[i+n,:] == Xs[j+n,:]
				push!(ii,j)
			end
		end
	end
	
	jj = setdiff(1:n,ii)

	ids = vec(readdlm("data_PEN/pen_$(date)_ids.csv",','))

	IDS = ids[jj]

	writedlm("data_PEN/pen_$(date)_final_ids.csv",IDS,',')
end


