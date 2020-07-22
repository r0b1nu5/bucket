using Dates

include("final_new_par.jl")

@info "====================================================================="
@info "$(now())"
@info "====================================================================="

ids = [
       "ntw3_1",
       "ntw3_2",
       "ntw3_3",
#       "ntw20_1",
#       "ntw20_2",
#       "ntw20_3",
#       "ntw20_4",
#       "ieee57_1",
#       "ieee57_2",
#       "pen_1",
#       "pen_2",
#       "pen_2",
#       "pen_3",
#       "pen_4",
#       "pen_5",
#       "pen_6",
#       "pen_7",
#       "pen_8"
      ]

ks_ntw3 = (1,50,1)
ks_ntw20 = (1,30,1)
ks_ieee57 = (1,30,1)
ks_pen = (1500,11000,1000)

files = Dict{String,String}(
			    "ntw3_1" => "data_marc/Xs1.csv",
			    "ntw3_2" => "data_marc/Xs2.csv",
			    "ntw3_3" => "data_marc/Xs3.csv",
			    "ntw20_1" => "data_melvyn/time_series_20_0.100000_0.010000_0.005000.csv",
			    "ntw20_2" => "data_melvyn/time_series_20_0.100000_0.010000_0.010000.csv",
			    "ntw20_3" => "data_melvyn/time_series_20_0.100000_0.100000_0.005000.csv",
			    "ntw20_4" => "data_melvyn/time_series_20_0.100000_0.100000_0.010000.csv",
			    "ieee57_1" => "data_melvyn/time_series_57_0.100000_0.010000_0.010000.csv",
			    "ieee57_2" => "data_melvyn/time_series_57_0.100000_0.010000_0.100000.csv",
			    "pen_1" => "data_PEN/Xs_2013-01-15_00.csv",
			    "pen_2" => "data_PEN/Xs_2013-03-10_04.csv",
			    "pen_3" => "data_PEN/Xs_2013-04-03_02.csv",
			    "pen_4" => "data_PEN/Xs_2013-04-03_03.csv",
			    "pen_5" => "data_PEN/Xs_2013-04-03_07.csv",
			    "pen_6" => "data_PEN/Xs_2013-07-30_01.csv",
			    "pen_7" => "data_PEN/Xs_2013-07-30_04.csv",
			    "pen_8" => "data_PEN/Xs_2013-07-30_09.csv"
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
			    "pen_1" => 1/30,
			    "pen_2" => 1/30,
			    "pen_3" => 1/30,
			    "pen_4" => 1/30,
			    "pen_5" => 1/30,
			    "pen_6" => 1/30,
			    "pen_7" => 1/30,
			    "pen_8" => 1/30
			    )

Ks = Dict{String,Tuple{Int64,Int64,Int64}}(
					   "ntw3_1" => ks_ntw3,
					   "ntw3_2" => ks_ntw3,
					   "ntw3_3" => ks_ntw3,
					   "ntw20_1" => ks_ntw20,
					   "ntw20_2" => ks_ntw20,
					   "ntw20_3" => ks_ntw20,
					   "ntw20_4" => ks_ntw20,
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

for id in ids
	Xs = readdlm(files[id],',')
	nn = size(Xs)[1]
	n = Int(nn/2)

	xxx = run_new_l0_par(id,Xs,taus[id],(1,n,1),Ks[id],false)
end


@info "====================================================================="
@info "$(now())"
@info "====================================================================="

