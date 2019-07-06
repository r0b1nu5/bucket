using PyPlot 

include("generate_time_series.jl")
include("locate_forced.jl")

ns = [5,]

dt = .1
eps = .05
Ts = Array{Int64,1}()	

for n in ns
	@info "============= n = $n =============="
	
	Xs = readdlm("data/C$(n)_$(dt).csv",',')
	
	L = 2*diagm(0 => ones(n)) - diagm(1 => ones(n-1)) - diagm(-1 => ones(n-1)) - diagm(n-1 => ones(1)) - diagm(1-n => ones(1))
	m = ones(n)
	d = ones(n)
	
	err = 1000.
	T = 200
	
	while err > eps
		@info "$(now()) -- error = $err"
		T += 200
		
		X = Xs[:,201:T]
		
		sol = locate_forced(X,m,d,dt)
		
		err = maximum(abs.(L - sol.L))
	end
	
	push!(Ts,T)
end

figure()
PyPlot.plot(ns,Ts,"-o")









