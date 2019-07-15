using LinearAlgebra, DelimitedFiles, Distributions

function generate_time_series(ntw::String, L::Array{Float64,2}, m::Array{Float64,1}, d::Array{Float64,1}, T::Int64, dt::Float64, sig::Array{Float64,1})
	n = size(L)[1]
	script_id = rand(1:1000)
	
	Mi = diagm(0 => 1 ./ m)
	D = diagm(0 => d)
	
	Ad = [zeros(n,n) diagm(0 => ones(n)); -Mi*L -Mi*D]
	A = Ad*dt + diagm(0 => ones(2*n))
	
	B = diagm(0 => sqrt(dt)*sig./m)
	
	X0 = zeros(2*n)
	X = Array{Float64,2}(undef,2*n,0)
	
	t = 0
	surcount = 0
	
	while t < T
		surcount += 1
		subcount = 0
		while subcount < 1000 && t < T
			t += 1
			subcount += 1
			
			xi = rand(Normal(0,1),n)
			
			X = [X (A*X0 + [zeros(n);B*xi])]
			X0 = X[:,end]
		end
		writedlm("data/temp_$(script_id)_$(surcount).csv,",X,',')
		X = Array{Float64,2}(undef,2*n,0)
	end
	
	Xf = Array{Float64,2}(undef,2*n,0)
	for i in 1:surcount
		Xf = [Xf readdlm("data/temp_$(script_id)_$(i).csv",',')]
		rm("data/temp_$(script_id)_$(i).csv")
	end
	
	writedlm("data/"*ntw*"_$(T)_$(dt).csv",Xf,',')
	
	return Xf
end


function generate_forced_time_series(ntw::String, L::Array{Float64,2}, m::Array{Float64,1}, d::Array{Float64,1}, forcing::Tuple{Array{Float64,1}, Array{Float64,1}, Array{Float64,1}}, T::Int64, dt::Float64, sig::Array{Float64,1})
	n = size(L)[1]
	script_id = rand(1:1000)
	
	Mi = diagm(0 => 1 ./ m)
	D = diagm(0 => d)
	
	Ad = [zeros(n,n) diagm(0 => ones(n)); -Mi*L -Mi*D]
	A = Ad*dt + diagm(0 => ones(2*n))
	
	B = diagm(0 => sqrt(dt)*sig./m)
	
	c,f,phi = forcing
	
	X0 = zeros(2*n)
	X = Array{Float64,2}(undef,2*n,0)
	
	t = 0

	str = "data/"*ntw*"_forced_$(maximum(abs.(f)))_$(T)_$(dt).csv"
 ##=
	surcount = 0
	
	while t < T
		@info "t/T = $(t)/$(T)"
		surcount += 1
		subcount = 0
		while subcount < 1000 && t < T
			t += 1
			subcount += 1
			
			xi = rand(Normal(0,1),n)
			
			X = [X (A*X0 + [zeros(n);B*xi] + [zeros(n);c.*cos.(2*pi*dt*t*f .+ phi)])]
			X0 = X[:,end]
		end
		writedlm("data/temp_$(script_id)_$(surcount).csv",X,',')
		X = Array{Float64,2}(undef,2*n,0)
	end

	Xf = Array{Float64,2}(undef,2*n,0)
	for i in 1:surcount
		X = readdlm("data/temp_$(script_id)_$(i).csv",',')
		Xf = [Xf X]
		rm("data/temp_$(script_id)_$(i).csv")
	end

	writedlm("data/"*ntw*"_forced_$(maximum(abs.(f)))_$(T)_$(dt).csv",Xf,',')
# =#
		
 #=	
	while t < T
		t += 1
		
		if t%1000 == 0
			@info "t/T = $t/$T"
		end
		
		xi = rand(Normal(0,1),n)
		
		X0 = A*X0 + [zeros(n);B*xi] + [zeros(n);c.*cos.(2*pi*dt*t*f .+ phi)]
		
		open(str, "a") do io
			writedlm(io, X0', ',')
		end
	end
# =#

end


function generate_forced_inertialess_time_series(ntw::String, L::Array{Float64,2}, forcing::Tuple{Array{Float64,1}, Array{Float64,1}, Array{Float64,1}}, T::Int64, dt::Float64, sig::Array{Float64,1})
	n = size(L)[1]
	script_id = rand(1:1000)
	
	B = diagm(0 => sqrt(dt)*sig)
	
	I = diagm(0 => ones(n))
	
	c,f,phi = forcing
	
	@info "c,f,phi = ($(c[45]),$(f[45]),$(phi[45])"
	
	X0 = zeros(n)
	X = Array{Float64,2}(undef,n,0)
	
	t = 0
	
	str = "data/"*ntw*"_forced_inertialess_$(maximum(abs.(f)))_$(T)_$(dt).csv"
	
	surcount = 0
	
	while t < T
		@info "t/T = $(t)/$(T)"
		surcount += 1
		subcount = 0
		while subcount < 1000 && t < T
			t += 1
			subcount += 1
			
			xi = rand(Normal(0,1),n)
			
			X = [X ((I - dt*L)*X0 + dt*B*xi + dt*c.*cos.(2*pi*f*dt*t .+ phi))]
			X0 = X[:,end]
		end
		writedlm("data/temp_$(script_id)_$(surcount).csv",X,',')
		X = Array{Float64,2}(undef,n,0)
	end
	
	Xf = Array{Float64,2}(undef,n,0)
	for i in 1:surcount
		X = readdlm("data/temp_$(script_id)_$(i).csv",',')
		Xf = [Xf X]
		rm("data/temp_$(script_id)_$(i).csv")
	end
	
	writedlm(str,Xf,',')

end




