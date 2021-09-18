using PyPlot, DelimitedFiles

to_plot = [("naspi_1",1),("naspi_2",1),("naspi_3",1),("naspi_4",1),("naspi_5",1),("naspi_6",1),("naspi_7",1),("naspi_8",1),("naspi_9",1)]
to_plot = [("naspi_10",1),("naspi_11",1),("naspi_12",1),("naspi_13",1)]

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
							 ("pen_8",2) => (9000,10000,50),
							 ("naspi_1",1) => (1,100,1),
							 ("naspi_2",1) => (1,100,1),
							 ("naspi_3",1) => (1,100,1),
							 ("naspi_4",1) => (1,100,1),
							 ("naspi_5",1) => (1,100,1),
							 ("naspi_6",1) => (1,100,1),
							 ("naspi_7",1) => (1,100,1),
							 ("naspi_8",1) => (1,100,1),
							 ("naspi_9",1) => (1,100,1),
							 ("naspi_10",1) => (1,100,1),
							 ("naspi_11",1) => (1,100,1),
							 ("naspi_12",1) => (1,100,1),
							 ("naspi_13",1) => (1,100,1),
							 ("mysterious_forcing_UK",1) => (1,50,1),
							 ("mysterious_forcing_57",1) => (1,50,1),
							 ("mysterious_forcing",1) => (1,1500,1),
							 ("mysterious_forcing",2) => (1,1500,1),
							 ("ntw20_multisine",1) => (10,200,10),
							 ("ntw20_saw",1) => (10,200,10),
							 ("ntw20_step",1) => (10,200,10),
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
			"pen_8" => 134,
			"naspi_1" => 58,
			"naspi_2" => 58,
			"naspi_3" => 58,
			"naspi_4" => 58,
			"naspi_5" => 58,
			"naspi_6" => 58,
			"naspi_7" => 58,
			"naspi_8" => 58,
			"naspi_9" => 58,
			"naspi_10" => 51,
			"naspi_11" => 58,
			"naspi_12" => 58,
			"naspi_13" => 58,
			"mysterious_forcing_UK" => 120,
			"mysterious_forcing_57" => 57,
			"mysterious_forcing" => 99,
			"ntw20_multisine" => 20,
			"ntw20_saw" => 20,
			"ntw20_step" => 20,
			)

node_ids = Dict{String,Array{Any,1}}(
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
				       "pen_8" => vec(Int.(readdlm("data_pen/pen_2013-07-30_09_ids.csv",',').-1)), 
				       "naspi_1" => vec(readdlm("data_naspi/naspi_ids_Case1.csv",',')),
				       "naspi_2" => vec(readdlm("data_naspi/naspi_ids_Case1.csv",',')),
				       "naspi_3" => vec(readdlm("data_naspi/naspi_ids_Case1.csv",',')),
				       "naspi_4" => vec(readdlm("data_naspi/naspi_ids_Case1.csv",',')),
				       "naspi_5" => vec(readdlm("data_naspi/naspi_ids_Case1.csv",',')),
				       "naspi_6" => vec(readdlm("data_naspi/naspi_ids_Case1.csv",',')),
				       "naspi_7" => vec(readdlm("data_naspi/naspi_ids_Case1.csv",',')),
				       "naspi_8" => vec(readdlm("data_naspi/naspi_ids_Case1.csv",',')),
				       "naspi_9" => vec(readdlm("data_naspi/naspi_ids_Case1.csv",',')),
				       "naspi_10" => vec(readdlm("data_naspi/naspi_ids_Case1.csv",','))[1:51],
				       "naspi_11" => vec(readdlm("data_naspi/naspi_ids_Case1.csv",',')),
				       "naspi_12" => vec(readdlm("data_naspi/naspi_ids_Case1.csv",',')),
				       "naspi_13" => vec(readdlm("data_naspi/naspi_ids_Case1.csv",',')),
				       "mysterious_forcing_UK" => Array(1:120),
				       "mysterious_forcing_57" => Array(1:57),
				       "mysterious_forcing" => Array(1:99),
				       "ntw20_multisine" => Array(1:20),
				       "ntw20_saw" => Array(1:20),
				       "ntw20_step" => Array(1:20),
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
		       "pen_8" => 1260.,
		       "naspi_1" => 56.6, 
		       "naspi_2" => 1778/30,
		       "naspi_3" => 1758/30,
		       "naspi_4" => 1898/30,
		       "naspi_5" => 1698/30,
		       "naspi_6" => 1778/30,
		       "naspi_7" => 1848/30,
		       "naspi_8" => 1848/30,
		       "naspi_9" => 1788/30,
		       "naspi_10" => 1798/30,
		       "naspi_11" => 1798/30,
		       "naspi_12" => 1898/30,
		       "naspi_13" => 1598/30,
		       "mysterious_forcing_UK" => 50000*.01,
		       "mysterious_forcing_57" => 200000*2e-3,
		       "mysterious_forcing" => 3000*.1,
		       "ntw20_multisine" => 100.,
		       "ntw20_saw" => 100.,
		       "ntw20_step" => 100.,
		       )



Ls = Dict{String,Array{Float64,2}}()
cmap = get_cmap("plasma")

for ntw_run in to_plot
	ntw,run = ntw_run
	if ntw_run == ("mysterious_forcing",2)
		ls = [1:84;86:99]
		ni = ls
	else
		ls = 1:ns[ntw]
		ni = node_ids[ntw]
	end
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

	node_id = node_ids[ntw][iii[1]]
	freq = round(ks[iii[2]]/T[ntw],digits=3)
	
	sort_nodes = sortslices([L[:,iii[2]] ni Array(1:length(ls))],dims=1)
	
	figure("[ℓ0] ntw: "*ntw*", run: $run",(8.,7.))

#=
	subplot2grid((1,7),(0,0),colspan=3)
	for i in 1:length(ls)
		j = Int64(sort_nodes[i,3])
		PyPlot.plot(ks/T[ntw],L[j,:],"-o",label="l = $(ls[j])",color=cmap((i-1)/(length(ls)-1)))
	end
	xlabel("freq")
	ylabel("objective")
	twiny()
	PyPlot.plot([ks[1],ks[end]],[Lmi,Lmi],"--k")
	xlabel("k")
=#

#	subplot2grid((1,7),(0,3),colspan=3)
	subplot2grid((1,6),(0,0),colspan=5)
	for i in 1:length(ls)
		j = Int(sort_nodes[i,3])
		PyPlot.plot(ks/T[ntw],L[j,:],"-o",label="l = $(ls[j])",color=cmap((i-1)/(length(ls)-1)))
	end
	xlabel("freq")
	ylabel("objective")
	twiny()
	PyPlot.plot([ks[1],ks[end]],[Lmi,Lmi],"--k")
	xlabel("k")
	title("Node id: $(node_id) \n frequency: $freq [Hz]")
	
#	subplot2grid((1,7),(0,6),colspan=1)
	subplot2grid((1,6),(0,5),colspan=1)
	xticks([])
	yticks([])
	twinx()
	yticks(-ls,sort_nodes[:,2])
	for i in 1:length(ls)
		PyPlot.plot([0,1],[-i,-i],color=cmap((i-1)/(length(ls)-1)),linewidth=2)
#		PyPlot.text(0,-i,sort_nodes[i,2])
	end
	axis([-.1,1.1,-maximum(ls)-1,-minimum(ls)+1])
end

					   


