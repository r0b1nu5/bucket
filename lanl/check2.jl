using PyPlot 

include("generate_time_series.jl")
include("locate_forced.jl")

ns = [5,10,20,40]

dt = .01
eps = .1
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
		
		S0 = zeros(2*n,2*n)
		S1 = zeros(2*n,2*n)
		for i in 1:T-200
			S0 += Xs[:,200+i]*transpose(Xs[:,200+i]) ./ (T-200)
			S1 += Xs[:,200+i+1]*transpose(Xs[:,200+i]) ./ (T-200)
		end
		
		Ah = S1*inv(S0)
		
		Lh = Ah[n+1:2*n,1:n]
				
		err = maximum(abs.(L - Lh))
	end
	
	push!(Ts,T)
end

figure()
PyPlot.plot(ns,Ts,"-o")









