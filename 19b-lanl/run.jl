using Dates

# Choose between l0 and l1 penalty.
approach = "l0"
#approach = "l1"

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
#	"uk_1"
#	"ebc_1",
       "ebc_2",
#       "ebc_3",
#       "ebc_4",
#       "ebc_5",
#       "ebc_6",
#       "ebc_7",
#       "ebc_8",
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
ks_ebc = Array(6500:1:6600)
#ks_ebc = Array(1500:1000:11000)
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
			    "ebc_1" => "data_ebc/Xs_2013-01-15_00.csv",
			    "ebc_2" => "data_ebc/Xs_2013-03-10_04.csv",
			    "ebc_3" => "data_ebc/Xs_2013-04-03_02.csv",
			    "ebc_4" => "data_ebc/Xs_2013-04-03_03.csv",
			    "ebc_5" => "data_ebc/Xs_2013-04-03_07.csv",
			    "ebc_6" => "data_ebc/Xs_2013-07-30_01.csv",
			    "ebc_7" => "data_ebc/Xs_2013-07-30_04.csv",
			    "ebc_8" => "data_ebc/Xs_2013-07-30_09.csv",
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
			    "ebc_1" => 1/30,
			    "ebc_2" => 1/30,
			    "ebc_3" => 1/30,
			    "ebc_4" => 1/30,
			    "ebc_5" => 1/30,
			    "ebc_6" => 1/30,
			    "ebc_7" => 1/30,
			    "ebc_8" => 1/30,
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
					   "ebc_1" => ks_ebc,
					   "ebc_2" => ks_ebc,
					   "ebc_3" => ks_ebc,
					   "ebc_4" => ks_ebc,
					   "ebc_5" => ks_ebc,
					   "ebc_6" => ks_ebc,
					   "ebc_7" => ks_ebc,
					   "ebc_8" => ks_ebc,
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
		global xxx = run_l0_par(id,Xs,taus[id],Array(1:1:n),Ks[id],false)
	elseif approach == "l1"
		global xxx = run_l1_par(id,Xs,taus[id],Ks[id],false)
	end
end


@info "====================================================================="
@info "$(now())"
@info "====================================================================="
