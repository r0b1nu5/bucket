using Distributed

@everywhere include("load_pen.jl")
@everywhere include("final.jl")

@everywhere function loc_pen(date::String)
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
	
	Xs = Xs[[jj;jj.+n],:]

	dt = 1/30
	
	Lh,dh,ah,fh,ph = run_location_large_ntw(Xs, dt, 5)

	writedlm("data_PEN/pen_"*date*"_Lh.csv",Lh,',')
	writedlm("data_PEN/pen_"*date*"_dh.csv",dh,',')
	writedlm("data_PEN/pen_"*date*"_ah.csv",ah,',')
	writedlm("data_PEN/pen_"*date*"_fh.csv",fh,',')
	writedlm("data_PEN/pen_"*date*"_ph.csv",ph,',')
end

pmap(loc_pen,dates)


