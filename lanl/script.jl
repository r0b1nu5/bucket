using Dates 

include("generate_time_series.jl")

dt = .25
T = 50000

for n in [5,10,20,40,80]
	@info "n = $(n)"
	
	ntw = "C$(n)"
	
	L = 2*diagm(0 => ones(n)) - diagm(1 => ones(n-1)) - diagm(-1 => ones(n-1)) - diagm(n-1 => ones(1)) - diagm(1-n => ones(1))
	m = ones(n)
	d = ones(n)
	
	X = generate_time_series(ntw,L,m,d,T,dt,ones(n))
end





