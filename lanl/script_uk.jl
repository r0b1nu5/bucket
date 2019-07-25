using Dates 

include("generate_time_series.jl")


for i in [1, 3, 8, 9, 13, 14, 15, 17]
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
	
	generate_forced_time_series("uk$(i)",L,m,d,(a,f,p),100000,.1,ones(120))
	
	@info "$(now()) -- $i terminated."
end


