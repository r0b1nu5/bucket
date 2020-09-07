using PyPlot, DelimitedFiles

ids = [
       "ntw3_1",
       "ntw3_2",
       "ntw3_3",
       "ntw20_1",
       "ntw20_2",
       "ntw20_3",
       "ntw20_4",
       "ntw20_5",
       "ntw20_6",
       "ieee57_1",
       "ieee57_2",
       "pen_1",
       "pen_2",
       "pen_3",
       "pen_4",
       "pen_5",
       "pen_6",
       "pen_7",
       "pen_8"
      ]

ks_ntw3 = (1,50,1)
ks_ntw20 = (1,30,1)
ks_ieee57 = (1,30,1)
#ks_pen = (1500,11000,1000)	## 1st range for all pen
ks_pen = (5000,6000,50)	## 2nd range for pen_2 and pen_5
#ks_pen = (8500,10000,50)	## 2nd range for pen_3
#ks_pen = (6000,7000,50)		## 2nd range for pen_4
#ks_pen = (9000,10000,50)	## 2nd range for pen_8

Ks = Dict{String,Tuple{Int64,Int64,Int64}}(
					   "ntw3_1" => ks_ntw3,
					   "ntw3_2" => ks_ntw3,
					   "ntw3_3" => ks_ntw3,
					   "ntw20_1" => ks_ntw20,
					   "ntw20_2" => ks_ntw20,
					   "ntw20_3" => ks_ntw20,
					   "ntw20_4" => ks_ntw20,
					   "ntw20_5" => ks_ntw20,
					   "ntw20_6" => ks_ntw20,
					   "ieee57_1" => ks_ieee57,
					   "ieee57_2" => ks_ieee57,
					   "pen_1" => ks_pen,
					   "pen_2" => ks_pen,
					   "pen_3" => ks_pen,
					   "pen_4" => ks_pen,
					   "pen_5" => ks_pen,
					   "pen_6" => ks_pen,
					   "pen_7" => ks_pen,
					   "pen_8" => ks_pen
					   )

ns = Dict{String,Int64}(
			"ntw3_1" => 3,
			"ntw3_2" => 3,
			"ntw3_3" => 3,
			"ntw20_1" => 20,
			"ntw20_2" => 20,
			"ntw20_3" => 20,
			"ntw20_4" => 20,
			"ntw20_5" => 20,
			"ntw20_6" => 20,
			"ieee57_1" => 57,
			"ieee57_2" => 57,
			"pen_1" => 115,
			"pen_2" => 129,
			"pen_3" => 130,
			"pen_4" => 129,
			"pen_5" => 130,
			"pen_6" => 137,
			"pen_7" => 136,
			"pen_8" => 134
			)


Ls = Dict{String,Array{Float64,2}}()
#for ntw in ids
for ntw in ["pen_2",]
	ls = 1:ns[ntw]
	ks = Ks[ntw][1]:Ks[ntw][3]:Ks[ntw][2]
	L = zeros(length(ls),length(ks))
	for i in 1:length(ls)
		for j in 1:length(ks)
			L[i,j] = readdlm("data/"*ntw*"_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
		end
	end
	Ls[ntw] = L

	figure(ntw)
	for i in 1:length(ls)
		PyPlot.plot(ks,L[i,:],"-o",label="l = $(ls[i])")
	end
	xlabel("k")
	ylabel("objective")
end

					   


