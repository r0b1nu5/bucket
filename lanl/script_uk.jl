using Dates 

include("generate_time_series.jl")

 #= Inhomogeneities

for i in 11:20
	@info "$(now()) -- $i starts..."
	
	L,m,d = generate_uk_ntw((.5,2.),(.5,2.),(.5,2.))
	
	writedlm("data/uk$(i)_L.csv",L,',')
	writedlm("data/uk$(i)_m.csv",m,',')
	writedlm("data/uk$(i)_d.csv",d,',')
	
	a = zeros(120)
	a[45] = .2
	f = zeros(120)
	f[45] = .1
	p = zeros(120)
	p[45] = pi/6
	
	generate_forced_time_series("uk$(i)",L,m,d,(a,f,p),100000,.01,ones(120))
	
	@info "$(now()) -- $i terminated."
end
# =#

# #= dt

include("load_uk.jl")
dts = [.1,.05,.01,.005,.002,.001,.0005]

a = zeros(120)
a[45] = .2
f = zeros(120)
f[45] = .1
p = zeros(120)
p[45] = pi/6

for dt in dts
	@info "$(now()) -- δt = $(dt)"

	generate_forced_time_series("uk",L,m,d,(a,f,p),100000,dt,ones(120))
	
	@info "$(now()) -- δt = $(dt) terminated."
end
# =#

# #= T

include("load_uk.jl")
Ts = [1000, 5000, 10000, 20000, 50000, 75000, 100000, 200000]

a = zeros(120)
a[45] = .2
f = zeros(120)
f[45] = .1
p = zeros(120)
p[45] = pi/6

for T in Ts
	@info "$(now()) -- T = $(T)"
	
	generate_forced_time_series("uk",L,m,d,(a,f,p),T,.01,ones(120))

	@info "$(now()) -- T = $(T) terminated."
end
# =#

# #= dt*T

Tdts = [(10000,.1), (20000,.05), (50000,.02), (100000,.01), (200000,.005)]

a = zeros(120)
a[45] = .2
f = zeros(120)
f[45] = .1
p = zeros(120)
p[45] = pi/6

for (T,dt) in Tdts
	@info "$(now()) -- T = $(T), δt = $(dt)"

	generate_forced_time_series("uk",L,m,d,(a,f,p),T,dt,ones(120))

	@info "$(now()) -- T = $(T), δt = $(dt) terminated."
end
# =#











