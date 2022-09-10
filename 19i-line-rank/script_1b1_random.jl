include("1b1.jl")

 #=
ntw = "uk"
#Ps = [.5,.55,.6,.65,.7,.75,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.72,1.725]
Ps = [.65,.7,.75,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.72,1.725]
m = .2*ones(120)
d = .1*ones(120)
# =#
 #=
ntw = "ieee57"
#Ps = [.01,.02,.03,.04,.05,.06,.07,.08,.09,.1,.2,.3,.4,.5,.6,.65,.66,.67,.68,.684]
Ps = [.68,.684]
m = .2*ones(57)
d = .1*ones(57)
# =#
 #=
ntw = "ieee118"
#Ps = [.01,.015,.02,.025,.03,.035,.04,.041,.0412]
Ps = [.03,.035,.04,.041,.0412]
m = .2*ones(118)
d = .1*ones(118)
# =#
# #=
ntw = "ieee300"
Ps = [.002,.004,.006,.008,.01,.012,.014,.016,.018,.022,.024,.05]
m = 2*ones(300)
d = 1*ones(300)
# =#

meas = "Omega"

for P0 in Ps
	for i in 1:20
		ranks,rmvd,cuts = rmv_1b1(ntw,"random",meas,P0,m,d,)
		writedlm("data/random/ranks_1b1_"*ntw*"_$(P0)_random_$i.csv",ranks,',')
		writedlm("data/random/rmvd_1b1_"*ntw*"_$(P0)_random_$i.csv",rmvd,',')
		writedlm("data/random/cuts_1b1_"*ntw*"_$(P0)_random_$i.csv",cuts,',')
	end
end

