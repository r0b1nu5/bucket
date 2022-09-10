using PyPlot, DelimitedFiles

to_plot = [("naspi_1",1),("naspi_2",1),("naspi_3",1),("naspi_4",1),("naspi_5",1),("naspi_6",1),("naspi_7",1),("naspi_8",1),("naspi_9",1),("naspi_10",1),("naspi_11",1),("naspi_12",1),("naspi_13",1)]

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
							 ("ebc_1",1) => (1500,11000,1000),
							 ("ebc_2",1) => (1500,11000,1000),
							 ("ebc_2",2) => (5000,6000,50),
							 ("ebc_2",3) => (5750,5850,1),
							 ("ebc_3",1) => (1500,11000,1000),
							 ("ebc_3",2) => (8500,10000,50),
							 ("ebc_3",3) => (9300,9600,1),
							 ("ebc_4",1) => (1500,11000,1000),
							 ("ebc_4",2) => (6000,7000,50),
							 ("ebc_4",3) => (6500,6600,1),
							 ("ebc_5",1) => (1500,11000,1000),
							 ("ebc_5",2) => (5000,6500,50),
							 ("ebc_6",1) => (1500,11000,1000),
							 ("ebc_7",1) => (1500,11000,1000),
							 ("ebc_8",1) => (1500,11000,1000),
							 ("ebc_8",2) => (9000,10000,50),
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
							 ("naspi_a1",1) => (1,100,1),
							 ("naspi_a2",1) => (1,100,1),
							 ("naspi_a3",1) => (1,100,1),
							 ("mysterious_forcing_UK",1) => (1,50,1),
							 ("mysterious_forcing_57",1) => (1,50,1),
							 ("mysterious_forcing",1) => (1,1500,1),
							 ("mysterious_forcing_ebc2",1) => (5000,6000,50),
							 ("mysterious_forcing_ebc3",1) => (8500,10000,50),
							 ("mysterious_forcing_ebc4",1) => (6000,7000,50),
							 ("mysterious_forcing_ebc5",1) => (5000,6500,50),
							 ("mysterious_forcing_ebc8",1) => (9000,10000,50),
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
			"ebc_1" => 115,
			"ebc_2" => 129,
			"ebc_3" => 130,
			"ebc_4" => 129,
			"ebc_5" => 130,
			"ebc_6" => 137,
			"ebc_7" => 136,
			"ebc_8" => 134,
			"naspi_1" => 58,
			"naspi_2" => 58,
			"naspi_3" => 58,
			"naspi_4" => 58,
			"naspi_5" => 58,
			"naspi_6" => 58,
			"naspi_7" => 58,
			"naspi_8" => 58,
			"naspi_9" => 58,
			"naspi_10" => 58,
			"naspi_11" => 58,
			"naspi_12" => 58,
			"naspi_13" => 58,
			"naspi_a1" => 58,
			"naspi_a2" => 58,
			"naspi_a3" => 58,
			"mysterious_forcing_UK" => 120,
			"mysterious_forcing_57" => 57,
			"mysterious_forcing" => 99,
			"mysterious_forcing_ebc2" => 129,
			"mysterious_forcing_ebc3" => 130,
			"mysterious_forcing_ebc4" => 129,
			"mysterious_forcing_ebc5" => 130,
			"mysterious_forcing_ebc8" => 134,
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
				       "ebc_1" => vec(Int.(readdlm("data_ebc/ebc_2013-01-15_00_ids.csv",',').-1)),
				       "ebc_2" => vec(Int.(readdlm("data_ebc/ebc_2013-03-10_04_ids.csv",',').-1)),
				       "mysterious_forcing_ebc2" => vec(Int.(readdlm("data_ebc/ebc_2013-03-10_04_ids.csv",',').-1)),
				       "ebc_3" => vec(Int.(readdlm("data_ebc/ebc_2013-04-03_02_ids.csv",',').-1)),
				       "mysterious_forcing_ebc3" => vec(Int.(readdlm("data_ebc/ebc_2013-04-03_02_ids.csv",',').-1)),
				       "ebc_4" => vec(Int.(readdlm("data_ebc/ebc_2013-04-03_03_ids.csv",',').-1)),
				       "mysterious_forcing_ebc4" => vec(Int.(readdlm("data_ebc/ebc_2013-04-03_03_ids.csv",',').-1)),
				       "ebc_5" => vec(Int.(readdlm("data_ebc/ebc_2013-04-03_07_ids.csv",',').-1)),
				       "mysterious_forcing_ebc5" => vec(Int.(readdlm("data_ebc/ebc_2013-04-03_07_ids.csv",',').-1)),
				       "ebc_6" => vec(Int.(readdlm("data_ebc/ebc_2013-07-30_01_ids.csv",',').-1)),
				       "ebc_7" => vec(Int.(readdlm("data_ebc/ebc_2013-07-30_04_ids.csv",',').-1)),
				       "ebc_8" => vec(Int.(readdlm("data_ebc/ebc_2013-07-30_09_ids.csv",',').-1)),
				       "mysterious_forcing_ebc8" => vec(Int.(readdlm("data_ebc/ebc_2013-07-30_09_ids.csv",',').-1)),
				       "naspi_1" => vec(readdlm("data_naspi/naspi_ids_Case1.csv",',')),
				       "naspi_2" => vec(readdlm("data_naspi/naspi_ids_Case2.csv",',')),
				       "naspi_3" => vec(readdlm("data_naspi/naspi_ids_Case3.csv",',')),
				       "naspi_4" => vec(readdlm("data_naspi/naspi_ids_Case4.csv",',')),
				       "naspi_5" => vec(readdlm("data_naspi/naspi_ids_Case5.csv",',')),
				       "naspi_6" => vec(readdlm("data_naspi/naspi_ids_Case6.csv",',')),
				       "naspi_7" => vec(readdlm("data_naspi/naspi_ids_Case7.csv",',')),
				       "naspi_8" => vec(readdlm("data_naspi/naspi_ids_Case8.csv",',')),
				       "naspi_9" => vec(readdlm("data_naspi/naspi_ids_Case9.csv",',')),
				       "naspi_10" => vec(readdlm("data_naspi/naspi_ids_Case10.csv",',')),
				       "naspi_11" => vec(readdlm("data_naspi/naspi_ids_Case11.csv",',')),
				       "naspi_12" => vec(readdlm("data_naspi/naspi_ids_Case12.csv",',')),
				       "naspi_13" => vec(readdlm("data_naspi/naspi_ids_Case13.csv",',')),
				       "naspi_a1" => vec(readdlm("data_naspi/naspi_ids_A1.csv",',')),
				       "naspi_a2" => vec(readdlm("data_naspi/naspi_ids_A2.csv",',')),
				       "naspi_a3" => vec(readdlm("data_naspi/naspi_ids_A3.csv",',')),
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
		       "ebc_1" => 1565.6,
		       "ebc_2" => 1260.,
		       "ebc_3" => 660.,
		       "ebc_4" => 600.,
		       "ebc_5" => 1260.,
		       "ebc_6" => 660.,
		       "ebc_7" => 1260.,
		       "ebc_8" => 1260.,
		       "naspi_1" => 1698/30,
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
		       "naspi_a1" => 1598/30,
		       "naspi_a2" => 1598/30,
		       "naspi_a3" => 1598/30,
		       "mysterious_forcing_UK" => 50000*.01,
		       "mysterious_forcing_57" => 200000*2e-3,
		       "mysterious_forcing" => 3000*.1,
		       "mysterious_forcing_ebc2" => 37800/30,
		       "mysterious_forcing_ebc3" => 19800/30,
		       "mysterious_forcing_ebc4" => 18000/30,
		       "mysterious_forcing_ebc5" => 37800/30,
		       "mysterious_forcing_ebc8" => 37800/30,
		       "ntw20_multisine" => 100.,
		       "ntw20_saw" => 100.,
		       "ntw20_step" => 100.,
		       )



Ls = Dict{String,Array{Float64,1}}()
γs = Dict{String,Array{Array{Float64,1},1}}()
cmap = get_cmap("plasma")

for ntw_run in to_plot
	ntw,run = ntw_run
	Ks = kss[ntw_run]
	ks = Ks[1]:Ks[3]:Ks[2]
	L = Array{Float64,1}()
	γ = Array{Array{Float64,1},1}()
	for j in 1:length(ks)
		push!(L,readdlm("data/"*ntw*"_l1_$(ks[j])_obj.csv",',')[1])
		push!(γ,vec(readdlm("data/"*ntw*"_l1_$(ks[j])_g.csv",',')))
	end
	Ls[ntw] = L
	γs[ntw] = γ

	Lmi,iii = findmin(L)
	Lma = maximum(L)
	nma,jjj = findmax(γ[iii])

	node_id = node_ids[ntw][jjj]
	freq = round(ks[iii]/T[ntw],digits=3)
	
	figure("[ℓ1] ntw: "*ntw*", run: $run",(10.,4.))

	subplot2grid((1,21),(0,0),colspan=9)
	PyPlot.plot(ks/T[ntw],L,"-o",color=cmap(.4))
	xlabel("freq")
	ylabel("obj")
	twiny()
	PyPlot.plot([ks[1],ks[end]],[Lmi,Lmi],"--k")
	xlabel("k")

	subplot2grid((1,21),(0,11),colspan=9)
#	for i in 1:length(ks)
#		α = cmap(1 - (L[i] - Lma)/(Lmi - Lma))
#		PyPlot.plot(1:length(γ[i]),γ[i],color=α)
#	end
#	PyPlot.plot(jjj,γ[iii][jjj],"or",markersize=10)
#	PyPlot.plot(jjj,γ[iii][jjj],"ow",markersize=7)
	PyPlot.plot(1:length(γ[iii]),γ[iii],"-o",color=cmap(0.),label="Min objective, k = $(ks[iii])")
	legend()
	xlabel("node idx")
	ylabel("amplitude")
	title("Node id: $(node_id) \n frequency: $freq [Hz]")

#=
	subplot2grid((1,20),(0,19),colspan=1)
	yticks([])
	twinx()
	for ll in LinRange(Lmi,Lma,1000)
		α = cmap(1 - (ll - Lma)/(Lmi - Lma))
		PyPlot.plot([0,1],[ll,ll],color=α)
	end
	axis([0,1,Lmi,Lma])
	xticks([])
	ylabel("obj")
=#
end

					   


