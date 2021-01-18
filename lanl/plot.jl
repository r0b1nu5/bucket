using PyPlot, DelimitedFiles

to_plot = [("ntw3_1",1),("ntw3_2",1),("ntw3_3",1)]

kss = Dict{Tuple{String,Int64},Tuple{Int64,Int64,Int64}}(
							 ("ntw3_1",1) => (1,50,1),
							 ("ntw3_2",1) => (1,50,1),
							 ("ntw3_3",1) => (1,50,1),
							 ("ntw20_1",1) => (1,30,1),
							 ("ntw20_2",1) => (1,30,1),
							 ("ntw20_3",1) => (1,30,1),
							 ("ntw20_4",1) => (1,30,1),
							 ("ntw20_5",1) => (1,30,1),
							 ("ntw20_6",1) => (1,30,1),
							 ("ieee57_1",1) => (1,30,1),
							 ("ieee57_2",1) => (1,30,1),
							 ("pen_1",1) => (1500,11000,1000),
							 ("pen_2",1) => (1500,11000,1000),
							 ("pen_2",2) => (5000,6000,50),
							 ("pen_2",3) => (5750,5850,1),
							 ("pen_3",1) => (1500,11000,1000),
							 ("pen_3",2) => (8500,10000,50),
							 ("pen_3",3) => (9300,9600,1),
							 ("pen_4",1) => (1500,11000,1000),
							 ("pen_4",2) => (6000,7000,50),
							 ("pen_4",3) => (6500,6600,1),
							 ("pen_5",1) => (1500,11000,1000),
							 ("pen_5",2) => (5000,6500,50),
							 ("pen_6",1) => (1500,11000,1000),
							 ("pen_7",1) => (1500,11000,1000),
							 ("pen_8",1) => (1500,11000,1000),
							 ("pen_8",2) => (9000,10000,50)
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

node_ids = Dict{String,Array{Int64,1}}(
				       "ntw3_1" => [1,2,3],
				       "ntw3_2" => [1,2,3],
				       "ntw3_3" => [1,2,3],
				       "ntw20_1" => Array(1:20),
				       "ntw20_2" => Array(1:20),
				       "ntw20_3" => Array(1:20),
				       "ntw20_4" => Array(1:20),
				       "ntw20_5" => Array(1:20),
				       "ntw20_6" => Array(1:20),
				       "ieee57_1" => Array(1:57),
				       "ieee57_2" => Array(1:57),
				       "pen_1" => vec(Int.(readdlm("data_pen/pen_2013-01-15_00_ids.csv",',').-1)),
				       "pen_2" => vec(Int.(readdlm("data_pen/pen_2013-03-10_04_ids.csv",',').-1)),
				       "pen_3" => vec(Int.(readdlm("data_pen/pen_2013-04-03_02_ids.csv",',').-1)),
				       "pen_4" => vec(Int.(readdlm("data_pen/pen_2013-04-03_03_ids.csv",',').-1)),
				       "pen_5" => vec(Int.(readdlm("data_pen/pen_2013-04-03_07_ids.csv",',').-1)),
				       "pen_6" => vec(Int.(readdlm("data_pen/pen_2013-07-30_01_ids.csv",',').-1)),
				       "pen_7" => vec(Int.(readdlm("data_pen/pen_2013-07-30_04_ids.csv",',').-1)),
				       "pen_8" => vec(Int.(readdlm("data_pen/pen_2013-07-30_09_ids.csv",',').-1))
				       )

T = Dict{String,Float64}(
		       "ntw3_1" => 50.,
		       "ntw3_2" => 50.,
		       "ntw3_3" => 50.,
		       "ntw20_1" => 100.,
		       "ntw20_2" => 100.,
		       "ntw20_3" => 100.,
		       "ntw20_4" => 100.,
		       "ntw20_5" => 100.,
		       "ntw20_6" => 10.,
		       "ieee57_1" => 2000.,
		       "ieee57_2" => 20000.,
		       "pen_1" => 1565.6,
		       "pen_2" => 1260.,
		       "pen_3" => 660.,
		       "pen_4" => 600.,
		       "pen_5" => 1260.,
		       "pen_6" => 660.,
		       "pen_7" => 1260.,
		       "pen_8" => 1260.
		       )



Ls = Dict{String,Array{Float64,2}}()
for ntw_run in to_plot
	ntw,run = ntw_run
	ls = 1:ns[ntw]
	Ks = kss[ntw_run]
	ks = Ks[1]:Ks[3]:Ks[2]
	L = zeros(length(ls),length(ks))
	for i in 1:length(ls)
		for j in 1:length(ks)
			L[i,j] = readdlm("data/"*ntw*"_l0_$(ls[i]).$(ks[j])_obj.csv",',')[1]
		end
	end
	Ls[ntw] = L

	Lmi,iii = findmin(L)

#	node_id = iii[1]
	node_id = node_ids[ntw][iii[1]]
#	freq = ks[iii[2]]
	freq = ks[iii[2]]/T[ntw]
	
	figure("ntw: "*ntw*", run: $run")
	for i in 1:length(ls)
		PyPlot.plot(ks/T[ntw],L[i,:],"-o",label="l = $(ls[i])")
	end
	xlabel("freq")
	ylabel("objective")
	twiny()
	PyPlot.plot([ks[1],ks[end]],[Lmi,Lmi],"--k")
	xlabel("k")
	title("Node id: $(node_id), frequency: $freq")
end

					   


