using PyPlot 

include("generate_time_series.jl")
include("locate_forced.jl")

ns = [10,20]

dt = .3
eps = .06
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
	
	S0 = zeros(2*n,2*n)
	S1 = zeros(2*n,2*n)
	
	while err > eps
		@info "$(now()) -- T = $T, error = $err"
		
		T += 200
		
		for i in 1:200
			S0 += Xs[:,T-200+i]*transpose(Xs[:,T-200+i])
			S1 += Xs[:,T-200+i+1]*transpose(Xs[:,T-200+i])
		end
		
		Ah = S1*inv(S0)
		Adh = (Ah - diagm(0 => ones(2*n)))./dt
		
		err = maximum(abs.(Ad - Adh))
	end
	
	push!(Ts,T)
end

figure(1)
PyPlot.plot(ns,Ts,"-o",label="δt = $dt")









