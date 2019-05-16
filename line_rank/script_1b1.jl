include("1b1.jl")

 ##=
ntw = "uk"
Ps = [.5,.55,.6,.65,.7,.75,1.,1.1,1.2,1.3,1.4,1.5,1.6,1.7,1.72,1.725]
m = .2*ones(120)
d = .1*ones(120)
# =#
 #=
ntw = "ieee57"
Ps = [.01,.02,.03,.04,.05,.06,.07,.08,.09]
m = .2*ones(57)
d = .1*ones(57)
# =#
 #=
ntw = "ieee118"
Ps = [.01,.015,.02,.025,.03,.035,.04,.041,.0412]
m = .2*ones(118)
d = .1*ones(118)
# =#
 #=
ntw = "ieee300"
Ps = [.002,.004,.006,.008,.01,.012,.014,.016,.018,.022,.024,.05]
m = 2*ones(300)
d = 1*ones(300)
# =#
 #=
ntw = "pegase"
Ps = [.0025,.005,.0075,.01,.0125,.015,.0175,.02]
m = 2000*ones(2869)
d = 1000*ones(2869)
# =#

 #meas = "Omega"
 #meas = "Womega"
 #meas = "dKf1"
 #meas = "dKf2"
 #meas = "load"
 #meas = "Omega+load"

for meas in ["Womega","dKf2"]
for P0 in Ps
	rmv_1b1_2(ntw,"init",meas,P0,m,d)
	rmv_1b1_2(ntw,"updated",meas,P0,m,d)
#	rmv_1b1(ntw,"tini",meas,P0,m,d)
#	rmv_1b1(ntw,"detadpu",meas,P0,m,d)
end
end

