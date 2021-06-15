using Dates

# Choose between l0 and l1 penalty.
#approach = "l0"
approach = "l1"

include("final_par.jl")

@info "====================================================================="
@info "$(now())"
@info "====================================================================="

ids = [
#       "ntw3_1",
#       "ntw3_2",
#       "ntw3_3",
#       "ntw20_1",
#       "ntw20_2",
#       "ntw20_3",
#       "ntw20_4",
#       "ieee57_1",
#       "ieee57_2",
	"uk_1"
#	"pen_1",
#       "pen_2",
#       "pen_3",
#       "pen_4",
#       "pen_5",
#       "pen_6",
#       "pen_7",
#       "pen_8",
#       "naspi_1",
#       "naspi_2",
#       "naspi_3",
#       "naspi_4",
#       "naspi_5",
#       "naspi_6",
#       "naspi_7",
#       "naspi_8",
#       "naspi_9",
#       "naspi_10",
#       "naspi_11",
#       "naspi_12",
#       "naspi_13",
#       "naspi_a1",
#       "naspi_a2",
#       "naspi_a3",
]

ks_ntw3 = Array(1:1:50)
ks_ntw20 = Array(1:1:30)
ks_ieee57 = Array(1:1:30)
ks_uk = Array(1:1:30)
ks_pen = Array(6500:1:6600)
#ks_pen = Array(1500:1000:11000)
ks_naspi = Array(1:1:100)

files = Dict{String,String}(
			    "ntw3_1" => "data_marc/Xs1.csv",
			    "ntw3_2" => "data_marc/Xs2.csv",
			    "ntw3_3" => "data_marc/Xs3.csv",
			    "ntw20_1" => "data_melvyn/time_series_20_0.100000_0.010000_0.005000.csv",
			    "ntw20_2" => "data_melvyn/time_series_20_0.100000_0.010000_0.010000.csv",
			    "ntw20_3" => "data_melvyn/time_series_20_0.100000_0.100000_0.005000.csv",
			    "ntw20_4" => "data_melvyn/time_series_20_0.100000_0.100000_0.010000.csv",
			    "uk_1" => "data_melvyn/uk_Xs.csv",
			    "ieee57_1" => "data_melvyn/time_series_57_0.100000_0.010000_0.010000.csv",
			    "ieee57_2" => "data_melvyn/time_series_57_0.100000_0.010000_0.100000.csv",
			    "pen_1" => "data_pen/Xs_2013-01-15_00.csv",
			    "pen_2" => "data_pen/Xs_2013-03-10_04.csv",
			    "pen_3" => "data_pen/Xs_2013-04-03_02.csv",
			    "pen_4" => "data_pen/Xs_2013-04-03_03.csv",
			    "pen_5" => "data_pen/Xs_2013-04-03_07.csv",
			    "pen_6" => "data_pen/Xs_2013-07-30_01.csv",
			    "pen_7" => "data_pen/Xs_2013-07-30_04.csv",
			    "pen_8" => "data_pen/Xs_2013-07-30_09.csv",
			    "naspi_1" => "data_naspi/Xs_Case1.csv",
			    "naspi_2" => "data_naspi/Xs_Case2.csv",
			    "naspi_3" => "data_naspi/Xs_Case3.csv",
			    "naspi_4" => "data_naspi/Xs_Case4.csv",
			    "naspi_5" => "data_naspi/Xs_Case5.csv",
			    "naspi_6" => "data_naspi/Xs_Case6.csv",
			    "naspi_7" => "data_naspi/Xs_Case7.csv",
			    "naspi_8" => "data_naspi/Xs_Case8.csv",
			    "naspi_9" => "data_naspi/Xs_Case9.csv",
			    "naspi_10" => "data_naspi/Xs_Case10.csv",
			    "naspi_11" => "data_naspi/Xs_Case11.csv",
			    "naspi_12" => "data_naspi/Xs_Case12.csv",
			    "naspi_13" => "data_naspi/Xs_Case13.csv",
			    "naspi_a1" => "data_naspi/Xs_A1.csv",
			    "naspi_a2" => "data_naspi/Xs_A2.csv",
			    "naspi_a3" => "data_naspi/Xs_A3.csv"
			    )

taus = Dict{String,Float64}(
			    "ntw3_1" => 1e-3,
			    "ntw3_2" => 1e-3,
			    "ntw3_3" => 1e-3,
			    "ntw20_1" => 1e-3,
			    "ntw20_2" => 1e-3,
			    "ntw20_3" => 1e-3,
			    "ntw20_4" => 1e-3,
			    "ieee57_1" => 1e-3,
			    "ieee57_2" => 1e-3,
			    "uk_1" => .01,
			    "pen_1" => 1/30,
			    "pen_2" => 1/30,
			    "pen_3" => 1/30,
			    "pen_4" => 1/30,
			    "pen_5" => 1/30,
			    "pen_6" => 1/30,
			    "pen_7" => 1/30,
			    "pen_8" => 1/30,
			    "naspi_1" => 1/30,
			    "naspi_2" => 1/30,
			    "naspi_3" => 1/30,
			    "naspi_4" => 1/30,
			    "naspi_5" => 1/30,
			    "naspi_6" => 1/30,
			    "naspi_7" => 1/30,
			    "naspi_8" => 1/30,
			    "naspi_9" => 1/30,
			    "naspi_10" => 1/30,
			    "naspi_11" => 1/30,
			    "naspi_12" => 1/30,
			    "naspi_13" => 1/30,
			    "naspi_a1" => 1/30,
			    "naspi_a2" => 1/30,
			    "naspi_a3" => 1/30
			    )

Ks = Dict{String,Array{Int64,1}}(
					   "ntw3_1" => ks_ntw3,
					   "ntw3_2" => ks_ntw3,
					   "ntw3_3" => ks_ntw3,
					   "ntw20_1" => ks_ntw20,
					   "ntw20_2" => ks_ntw20,
					   "ntw20_3" => ks_ntw20,
					   "ntw20_4" => ks_ntw20,
					   "ieee57_1" => ks_ieee57,
					   "ieee57_2" => ks_ieee57,
					   "uk_1" => ks_uk,
					   "pen_1" => ks_pen,
					   "pen_2" => ks_pen,
					   "pen_3" => ks_pen,
					   "pen_4" => ks_pen,
					   "pen_5" => ks_pen,
					   "pen_6" => ks_pen,
					   "pen_7" => ks_pen,
					   "pen_8" => ks_pen,
					   "naspi_1" => ks_naspi,
					   "naspi_2" => ks_naspi,
					   "naspi_3" => ks_naspi,
					   "naspi_4" => ks_naspi,
					   "naspi_5" => ks_naspi,
					   "naspi_6" => ks_naspi,
					   "naspi_7" => ks_naspi,
					   "naspi_8" => ks_naspi,
					   "naspi_9" => ks_naspi,
					   "naspi_10" => ks_naspi,
					   "naspi_11" => ks_naspi,
					   "naspi_12" => ks_naspi,
					   "naspi_13" => ks_naspi,
					   "naspi_a1" => ks_naspi,
					   "naspi_a2" => ks_naspi,
					   "naspi_a3" => ks_naspi
					   )

for id in ids
	Xs = readdlm(files[id],',')
	nn = size(Xs)[1]
	n = Int(nn/2)

	if approach == "l0"
		xxx = run_l0_par(id,Xs,taus[id],Array(1:1:n),Ks[id],false)
	elseif approach == "l1"
		xxx = run_l1_par(id,Xs,taus[id],Ks[id],false)
	end
end


@info "====================================================================="
@info "$(now())"
@info "====================================================================="
