using PyPlot, DelimitedFiles, RecipesBase, Shapefile, Colors

include("LoadTsatTxt.jl")
#include("data_naspi/osc_init.jl")
#include("data_naspi/naspi_files.jl")

function ebc_preprocess_data(ntw::String)
	n = 0
	ks = 0:0
	T = 0.
	file = "none"
	date = "none"
	xid = Int64[]
	if ntw == "ebc_2"
		n = 129
		ks = Array(5000:50:6000)
		T = 1260. 
		file = "ebc_2"
		date = "2013-03-10_04"
	elseif ntw =="ebc_22"
		n = 129
		ks = Array(5000:50:6000)
		T = 1260.
		file = "ebc_22"
		date = "2013-03-10_04"
	elseif ntw == "ebc_3"
		n = 130
		ks = Array(8500:50:10000)
		T = 660.
		file = "ebc_3"
		date = "2013-04-03_02"
		xid = [46,]
	elseif ntw == "ebc_4"
		n = 129
		ks = Array(6000:50:7000)
		T = 600.
		file = "ebc_4"
		date = "2013-04-03_03"
		xid = [46,] 
	elseif ntw == "ebc_5"
		n = 130
		ks = Array(5000:50:6500)
		T = 1260.
		file = "ebc_5"
		date = "2013-04-03_07"
		xid = [46,] 
	elseif ntw == "ebc_55"
		#n = 130
		#ks = Array(5000:50:6500)
		n = 129
		ks = Array(5000:50:6000)
		T = 1260.
		file = "ebc_55"
		date = "2013-04-03_07"
		xid = [46,]
	elseif ntw == "ebc_8"
		n = 134
		ks = Array(9000:50:10000)
		T = 1260.
		file = "ebc_8"
		date = "2013-07-30_09"
	else
		@info "Invalid entry."
	end
	
	ls = 1:n
	
	L0 = zeros(length(ls),length(ks))
	for i in 1:length(ls)
		for j in 1:length(ks)
			L0[i,j] = readdlm("data/ebc/"*file*"_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
		end
	end
	L0 = L0[setdiff(1:n,xid),:]
	nL0 = (L0 .- maximum(L0))./(maximum(L0) - minimum(L0))
	
	mi,pos = findmin(nL0)
	j = pos[1]
	k = pos[2]

	return n-length(xid), ks, T, file, date, nL0, j, k
end

function utk_preprocess_boundaries()
	states = Shapefile.Table("data_utk/states_shapes/tl_2021_us_state.shp")
	states = Shapefile.Table("data_utk/states_shapes/cb_2018_us_state_5m.shp")
	skip = ["California","Oregon","Washington","Arizona","Nevada","Utah","Idaho","Alaska","Hawaii","American Samoa","Guam","Commonwealth of the Northern Mariana Islands","Puerto Rico","United States Virgin Islands"]

	can = Shapefile.Table("data_utk/states_shapes/gpr_000a11a_e.shp")
	toplot = ["Quebec","Ontario","New Brunswick","Nova Scotia"]

	ss = empty(Shapefile.shapes(states))
	s_xy = Vector{Tuple{Vector{Float64},Vector{Float64}}}()
	red_i = Vector{Vector{Vector{Int64}}}()
	name = Vector{String}()
	for row in states
		if !(row.NAME in skip)
			push!(ss, Shapefile.shape(row))
			x = [Shapefile.shape(row).points[i].x for i in 1:length(Shapefile.shape(row).points)]
			y = [Shapefile.shape(row).points[i].y for i in 1:length(Shapefile.shape(row).points)]
			push!(s_xy,(x,y))
			push!(name,row.NAME)
			if row.NAME == "Ohio"
				push!(red_i,[32:length(x),])
			elseif row.NAME == "Alabama"
				push!(red_i,[40:length(x),])
			elseif row.NAME == "Wisconsin"
				push!(red_i,[183:length(x),])
			elseif row.NAME == "Mississippi"
				push!(red_i,[60:length(x),])
			elseif row.NAME == "North Carolina"
				push!(red_i,[128:length(x),])
			elseif row.NAME == "New York"
				push!(red_i,[114:length(x),])
			elseif row.NAME == "Rhode Island"
				push!(red_i,[183:length(x),])
			elseif row.NAME == "Virginia"
				push!(red_i,[205:length(x),])
			elseif row.NAME == "Louisiana"
				push!(red_i,[141:length(x),])
			elseif row.NAME == "Massachusetts"
				push!(red_i,[126:length(x),])
			elseif row.NAME == "Florida"
				push!(red_i,[149:length(x),])
			elseif row.NAME == "Maryland"
				push!(red_i,[82:length(x),])
			elseif row.NAME == "Maine"
				push!(red_i,[182:length(x),])
			elseif row.NAME == "Michigan"
				push!(red_i,[(218:1268),(1372:length(x))])
			else
				push!(red_i,[1:length(x),])
			end
		end
	end

	for row in can
		if (row.PRENAME in toplot)
			push!(ss, Shapefile.shape(row))
			x = [Shapefile.shape(row).points[i].x for i in 1:length(Shapefile.shape(row).points)]
			y = [Shapefile.shape(row).points[i].y for i in 1:length(Shapefile.shape(row).points)]
			push!(s_xy,(x,y))
			push!(name,row.PRENAME)
			if row.PRENAME == "Quebec"
				push!(red_i,[173:length(x),])
			else
				push!(red_i,[1:length(x),])
			end
		end
	end
	
	l_xy = Vector{Tuple{Vector{Float64},Vector{Float64}}}()
	lakes = ["Erie","Huron","Michigan","Ontario","StClair","Superior"]
	lids = [Vector(1:349),
		Vector(1:431),
		Vector(1:1345),
		Vector(1:659),
		Vector(1:1208),
		Vector(1:590),
		Vector(1:273),
		Vector(1:83),
		Vector(1:75),
		Vector(1:976),
		Vector(1:717)]
	for lake in lakes
		lak = Shapefile.Table("data_utk/states_shapes/hydro_p_Lake"*lake*".shp")
		for colak in lak
			sha = Shapefile.shape(colak)
			x = Vector{Float64}()
			y = Vector{Float64}()
			for p in sha.points
				push!(x,p.x)
				push!(y,p.y)
			end
			push!(l_xy,(x,y))
		end
	end

	return s_xy, red_i, name, l_xy, lids, lakes
end

function rm_nan(Xs::Array{Float64,2})
	n,T = size(Xs)

	Xt = Array{Float64,2}(undef,0,T)

	for i in 1:n
		x = Xs[i,:]
		xnan = isnan.(x)

		b_nan = Array{Int64,1}()
		e_nan = Array{Int64,1}()
		now_nan = false

		if isnan(x[1])
			now_nan = true
			push!(b_nan,1)
		end

		for j in 2:length(x)
			if now_nan && xnan[j-1:j] == [1,0]
				push!(e_nan,j-1)
				now_nan = false
			elseif !now_nan && xnan[j-1:j] == [0,1]
				push!(b_nan,j)
				now_nan = true
			end
		end

		if now_nan
			push!(e_nan,T)
		end

		for j in 1:length(b_nan)
			b = b_nan[j]
			e = e_nan[j]

			if b == 1
				x[b:e] = x[e+1]*ones(e-b+1)
			elseif e == T
				x[b:e] = x[b-1]*ones(e-b+1)
			else
				x[b:e] = (x[e+1]-x[b-1])/((e+1)-(b-1))*((b:e) .- (b-1)) .+ x[b-1]
			end
		end

		Xt = [Xt;x']
	end

	return Xt
end

# Estimates the derivatives of the time series x, with step size h, using the "five point stencil" method.

function stencil_5p(x::Array{Float64,1}, h::Float64)
	T = length(x)

	dx = Array{Float64,1}()
	for i in 3:T-2
		push!(dx,(x[i-2] - 8*x[i-1] + 8*x[i+1] - x[i+2])/(12*h))
	end

	return x[3:T-2],dx
end

function stencil_5p(X::Array{Float64,2}, h::Float64)
	n,T = size(X)

	dX = Array{Float64,2}(undef,0,T-4)
	for i in 1:n
		x,dx = stencil_5p(X[i,:],h)
		dX = [dX;dx']
	end

	return X[:,3:T-2],dX
end


function prepare_naspi(case::String, file::String, init::Int64)
	t,Xtt,tit = LoadTsatTxt(file*case*"/BusVolAng.txt")
	Xt = rm_nan(Xtt)

	τ = t[2]-t[1]
	
	X,dX = stencil_5p(Xt[:,init:end]*π/180,τ)

	Xs = [X;dX]

	writedlm("data_naspi/naspi_ids_"*case*".csv",tit,',')
	writedlm("data_naspi/Xs_"*case*".csv",Xs,',')

	return Xs
end
