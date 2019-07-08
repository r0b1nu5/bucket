using PyPlot 

include("generate_time_series.jl")
include("locate_forced.jl")

ns = [5,10,20,40,80]

dt = .2
eps = .1
Ts = Array{Int64,1}()	

for n in ns
	@info "============= n = $n =============="
	
	Xs = readdlm("data/C$(n)_50000_$(dt).csv",',')
	
	L = 2*diagm(0 => ones(n)) - diagm(1 => ones(n-1)) - diagm(-1 => ones(n-1)) - diagm(n-1 => ones(1)) - diagm(1-n => ones(1))
	m = ones(n)
	Mi = diagm(0 => 1 ./ m)
	d = ones(n)
	D = diagm(0 => d)
	Ad = [zeros(n,n) diagm(0 => ones(n));-Mi*L -Mi*D]
	
	err = 1000.
	T = 200
	
	while err > eps
		@info "$(now()) -- T = $T, error = $err"
	
		T += 200
		
		X = Xs[:,201:T]
		
		sol = locate_forced(X,m,d,dt)
		Adh = [zeros(n,n) diagm(0 => ones(n));-Mi*sol.L -Mi*D]
		
		err = maximum(abs.(Ad - Adh))
	end
	
	push!(Ts,T)
end

figure(1)
PyPlot.plot(ns,Ts,"--s",label="δt = $dt")









