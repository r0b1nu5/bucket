using Distributed

@everywhere include("load_pen.jl")
@everywhere include("final.jl")

@everywhere function loc_pen(date::String)
	Xs = load_pen(date)

	dt = 1/30

	Lh,dh,ah,fh,ph = run_location_large_ntw(Xs, dt, 5)

	writedlm("data_PEN/pen_"*date*"_Lh.csv",Lh,',')
	writedlm("data_PEN/pen_"*date*"_dh.csv",dh,',')
	writedlm("data_PEN/pen_"*date*"_ah.csv",ah,',')
	writedlm("data_PEN/pen_"*date*"_fh.csv",fh,',')
	writedlm("data_PEN/pen_"*date*"_ph.csv",ph,',')
end

pmap(loc_pen,dates)


