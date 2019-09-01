using SpecialFunctions

include("mle.jl")

# Our goodness-of-fit test follows the recommandation of Clauset09. Namely, using the fitted paramters for each law, we generate some synthetic data and then compare the real data to the synthetic data. Our p-value is the proportion of synthetic data whose Kolmogorov-Smirnov statistic is larger than the one of our real data. If p > 0.1, we assume that the distribution if good.


## ======================= Goodness-of-fit estimates ================================

function new_gof_pl(j::String, x::Array{Int64,2}, s0::Float64, C0::Float64, mi::Int64, n_sample::Int64=2500)
	KS0 = new_KS_pl(x,s0,C0)
	n_data = sum(x[2,:])
	
	KSs = Array{Float64,1}()
	for i in 1:n_sample
#		if (i%100 == 0) || i == n_sample
			@info "$(now()) -- "*j*", GoF power law: $i/$n_sample"
#		end
		z = new_rand_pl(rand(n_data), s0, C0, mi)
		s = new_mle_pl(z)
		C = 1/zeta(s,mi)
		KS = new_KS_pl(z,s,C)
		push!(KSs,KS)
	end
	
	p = sum(KSs .> KS0)/n_sample
	
	writedlm("analysis/"*j*"_KS0_pl.csv",KS0,',')
	writedlm("analysis/"*j*"_KSs_pl.csv",KSs,',')
	
	return p
end


function new_gof_plc(j::String, x::Array{Int64,2}, a0::Float64, l0::Float64, C0::Float64, mi::Int64, n_sample::Int64=2500)
	KS0 = new_KS_plc(x,a0,l0,C0)
	n_data = sum(x[2,:])
	
	KSs = Array{Float64,1}()
	for i in 1:n_sample
#		if (i%100 == 0) || i == n_sample
			@info "$(now()) -- "*j*", GoF power law with cutoff: $i/$n_sample"
#		end
		z = new_rand_plc(rand(n_data), a0, l0, C0, mi)
# We should do as commented, but due to computation time, we truncate the tail of our synthetic data. Thus the MLE of the parameters is quite bad (at least for small data set, whereas the synthetic data follow correctly our real data set and its fit. We then compute KS with respect to the estimated parameter a0 and s0. In this case it will give more accurate results.
#= 
		a,l = new_mle_plc(z)
		C = 1/(real(polylog(a,Complex(exp(-l)))) - sum((1:mi-1).^(-a).*exp.(-l*(1:mi-1))))
		KS = new_KS_plc(z,a,l,C)
=#
		KS = new_KS_plc(z,a0,l0,C0)
		push!(KSs,KS)
	end
	
	p = sum(KSs .> KS0)/n_sample
	
	writedlm("analysis/"*j*"_KS0_plc.csv",KS0,',')
	writedlm("analysis/"*j*"_KSs_plc.csv",KSs,',')
	
	return p
end


function new_gof_yule(j::String, x::Array{Int64,2}, a0::Float64, C0::Float64, mi::Int64, n_sample::Int64=2500)
	KS0 = new_KS_yule(x,a0,C0)
	n_data = sum(x[2,:])
	
	KSs = Array{Float64,1}()
	for i in 1:n_sample
#		if (i%100 == 0) || i == n_sample
			@info "$(now()) -- "*j*", GoF Yule law: $i/$n_sample"
#		end
		z = new_rand_yule(rand(n_data), a0, C0, mi)
		a = new_mle_yule(z,mi)
		C = 1/(1-(a-1)*sum(beta.(1:(mi-1),a)))
		KS = new_KS_yule(z,a,C)
		push!(KSs,KS)
	end

@info "$(maximum(KSs)), $KS0"
	p = sum(KSs .> KS0)/n_sample
	
	writedlm("analysis/"*j*"_KS0_yule.csv",KS0,',')
	writedlm("analysis/"*j*"_KSs_yule.csv",KSs,',')
	
	return p
end


function new_gof_exp(j::String, x::Array{Int64,2}, b0::Float64, C0::Float64, mi::Int64, n_sample::Int64=2500)
	KS0 = new_KS_exp(x,b0,C0)
	n_data = sum(x[2,:])
	KSs = Array{Float64,1}()
	for i in 1:n_sample
		if (i%100 == 0) || i == n_sample
			@info "$(now()) -- "*j*", GoF exponential distribution: $i/$n_sample"
		end
		z = new_rand_exp(rand(n_data), b0, C0, mi)
		b = new_mle_exp(z,mi)
		C = (1 - exp(-b))/exp(-b*mi)
		KS = new_KS_exp(z,b,C)
		push!(KSs,KS)
	end

	p = sum(KSs .> KS0)/n_sample
@info "$p"
	@info "$(maximum(KSs))"
	@info "$KS0"
	writedlm("analysis/"*j*"_KS0_exp.csv",KS0,',')
	writedlm("analysis/"*j*"_KSs_exp.csv",KSs,',')

	return p
end


function new_gof_poisson(j::String, x::Array{Int64,2}, m0::Float64, C0::Float64, mi::Int64, n_sample::Int64=2500)
	KS0 = new_KS_poisson(x,m0,C0)
	n_data = sum(x[2,:])

	KSs = Array{Float64,1}()
	for i in 1:n_sample
		if (i%100 == 0) || i == n_sample
			@info "$(now()) -- "*j*", GoF Poisson distribution: $i/$n_sample"
		end
		z = new_rand_poisson(rand(n_data), m0, C0, mi)
		m = new_mle_poisson(z,mi)
		C = 1/(exp(m) - sum((m.^(0:mi-1))./(factorial.(0:mi-1))))
		KS = new_KS_poisson(z,m,C)
		push!(KSs,KS)
	end

	p = sum(KSs .> KS0)/n_sample

	writedlm("analysis/"*j*"_KS0_poisson.csv",KS0,',')
	writedlm("analysis/"*j*"_KSs_poisson.csv",KSs,',')

	return p
end


## ============================ Kolmogorov-Smirnov measures ==================================

function new_KS_pl(x::Array{Int64,2}, s::Float64, C::Float64)
	mi = x[1,1]
	ma = x[1,end]
	n = sum(x[2,:])
	
	cdf = [C*mi^(-s),]
	ecdf = [x[2,1]/n,]
	count = 2
	for i in (mi+1):ma
		if i == x[1,count]
			count += 1
		end
		push!(cdf,cdf[end] + C*i^(-s))
		push!(ecdf,sum(x[2,1:count-1])/n)
	end
	
	KS = maximum(abs.(cdf - ecdf))
	
	return KS
end

function new_KS_plc(x::Array{Int64,2}, a::Float64, l::Float64, C::Float64)
	mi = x[1,1]
	ma = x[1,end]
	n = sum(x[2,:])
	
	cdf = [C*mi^(-a)*exp(-l*mi),]
	ecdf = [x[2,1]/n,]
	count = 2
	for i in (mi+1):ma
		if i == x[1,count]
			count += 1
		end
		push!(cdf,cdf[end] + C*i^(-a)*exp(-l*i))
		push!(ecdf,sum(x[2,1:count-1])/n)
	end
	
	KS = maximum(abs.(cdf - ecdf))
	
	return KS
end

function new_KS_yule(x::Array{Int64,2}, a::Float64, C::Float64)
	mi = x[1,1]
	ma = x[1,end]
	n = sum(x[2,:])
	
	cdf = [C*(a-1)*beta(mi,a),]
	ecdf = [x[2,1]/n,]
	count = 2
	for i in (mi+1):ma
		if i == x[1,count]
			count += 1
		end
		push!(cdf, cdf[end] + C*(a-1)*beta(i,a))
		push!(ecdf, sum(x[2,1:count-1])/n)
	end
	
	KS = maximum(abs.(cdf - ecdf))
	
	return KS
end

function new_KS_exp(x::Array{Int64,2}, b::Float64, C::Float64)
	mi = x[1,1]
	ma = x[1,end]
	n = sum(x[2,:])

	cdf = [C*exp(-b*mi),]
	ecdf = [x[2,1]/n,]
	count = 2
	for i in (mi+1):ma
		if i == x[1,count]
			count += 1
		end
		push!(cdf, cdf[end] + C*exp(-b*i))
		push!(ecdf, sum(x[2,1:count-1])/n)
	end

	KS = maximum(abs.(cdf - ecdf))

	return KS
end

function new_KS_poisson(x::Array{Int64,2}, m::Float64, C::Float64)
	mi = x[1,1]
	ma = x[1,end]
	n = sum(x[2,:])

	cdf = [C * m^mi/factorial(mi),]
	ecdf = [x[2,1]/n,]
	count = 2
	for i in (mi+1):ma
		if i == x[1,count]
			count += 1
		end
		push!(cdf, cdf[end] + C*m^i/factorial(Float64(i)))
		push!(ecdf, sum(x[2,1:count-1])/n)
	end

	KS = maximum(abs.(cdf - ecdf))

	return KS
end


## ================================== Synthetic variables generation ================================

# Returns the histogram, takes too long to generate the data. The error in the tail of the CDF is eps.

function new_rand_pl(yy::Array{Float64,1}, s::Float64, C::Float64, mi::Int64 = 1, eps::Float64 = 1e-4)
	x = 0.
	n = Int(mi - 1)
	
	y = sort(yy)./C
	ly = length(y)
		
	val_num = Array{Int64,2}(undef,2,0)
	
	id = 0
	
	while (ly - id)*(1 - C*x) > eps
		while y[id+1] > x 
			n += 1
			x += n^(-s)
		end
		di = sum(y[id+1:end] .< x)
		id += di
		val_num = [val_num [n,di]]
	end
	
	val_num = [val_num [n+1,(ly-id)]]
	
 #= # Takes far too much time
	z = Array{Int64,1}()
	for i in 1:length(ns)
		z = [z;i*ones(ns[i])]
	end
# =#
	
	return val_num
end

function new_rand_plc(yy::Array{Float64,1}, a::Float64, l::Float64, C::Float64, mi::Int64 = 1, eps::Float64 = 1e-4)
	x = 0.
	n = Int(mi - 1)
	
	y = sort(yy)./C
	ly = length(y)
	
	val_num = Array{Int64,2}(undef,2,0)
	
	id = 0
	
	while (ly - id)*(1 - C*x) > eps
		while y[id+1] > x
			n += 1
			x += n^(-a)*exp(-l*n)
		end
		di = sum(y[id+1:end] .< x)
		id += di
		val_num = [val_num [n,di]]
	end
	
	val_num = [val_num [n+1,(ly-id)]]
	
	return val_num
end


function new_rand_yule(yy::Array{Float64,1}, a::Float64, C::Float64, mi::Int64 = 1, eps::Float64 = 1e-4)
	x = 0.
	n = Int(mi - 1)
	
	y = sort(yy)./C
	ly = length(y)
	
	val_num = Array{Int64,2}(undef,2,0)
	
	id = 0
	
	while (ly - id)*(1 - C*x) > eps
		while y[id+1] > x
			n += 1
			x += (a-1)*beta(n,a)
		end
		di = sum(y[id+1:end] .< x)
		id += di
		val_num = [val_num [n,di]]
	end
	
	val_num = [val_num [n+1, (ly-id)]]
	
	return val_num
end


function new_rand_exp(yy::Array{Float64,1}, b::Float64, C::Float64, mi::Int64=1, eps::Float64=1e-4)
	x = 0.
	n = Int(mi - 1)

	y = sort(yy)./C
	ly = length(y)

	val_num = Array{Int64,2}(undef,2,0)

	id = 0

	while (ly - id)*(1 - C*x) > eps
		while y[id+1] > x
			n += 1
			x += exp(-b*n)
		end
		di = sum(y[id+1:end] .< x)
		id += di
		val_num = [val_num [n,di]]
	end

	val_num = [val_num [n+1, (ly-id)]]

	return val_num
end


function new_rand_poisson(yy::Array{Float64,1}, m::Float64, C::Float64, mi::Int64=1, eps::Float64=1e-4)
	x = 0.
	n = Int(mi - 1)

	y = sort(yy)./C
	ly = length(y)

	val_num = Array{Int64,2}(undef,2,0)

	id = 0
	di = 0

	while (ly - id)*(1 - C*x) > eps
		while y[id+1] > x
			n += 1
			x += m^n/factorial(Float64(n))
		end
		di = sum(y[id+1:end] .< x)
		id += di
		val_num = [val_num [n,di]]
	end

	val_num = [val_num [n+1,di]]

	return val_num
end


## Old functions ====================================================================
# ===================================================================================

function my_gof(x::Array{Float64,1}, n_bins::Int64, n_sample::Int64=2500)
	mi = minimum(x)
	ma = maximum(x)
	n_data = length(x)
	
	bins = exp(collect(linspace(log(mi-.5),log(ma+.5),nbins+1)))
	midbins = sqrt(bins[2:end].*bins[1:end-1])
	
	xx = Array{Float64,1}()
	for X in x
		push!(xx,maximum(midbins.*(bins[1:end-1] .<= X)))
	end
	
	s0 = ml_cautet(xx)
	C0 = 1/sum((mi:1e7).^(-s0))
	KS0 = my_KS(xx,s0,C0,bins,midbins)
	
	KSs = Array{Float64,1}()
	for i in 1:n_sample
		@info(" $i/$n_sample")
		z = Array{Float64,1}()
		for j in 1:n_data
			push!(z,pl(rand(),s0,C0,mi))
		end
		new_bins = copy(bins)
		while maximum(new_bins) < maximum(z)
			push!(new_bins,new_bins[end]^2/new_bins[end-1])
		end
		new_midbins = sqrt(new_bins[2:end].*new_bins[1:end-1])
		
		zz = Array{Float64,1}()
		for Z in z
			push!(zz,maximum(new_midbins.*(new_bins[1:end-1] .<= Z)))
		end
		
		s = ml_cautet(zz)
		C = 1/sum((mi:1e7).^(-s))
		KS = my_KS(z,s,C,new_bins,new_midbins)
		push!(KSs,KS)
	end
	
	p = sum(KSs .> KS0)/n_sample
	
	return p
end

	
function gof(x::Array{Float64,1}, n_sample::Int64 = 2500)
	mi = minimum(x)
	ma = maximum(x)
	s0 = ml_zipf(x)
	C0 = 1/sum((mi:1e7).^(-s0))
	KS0, dtail_max, dtail_min, dtail_2 = KS_pl(x,s0,C0)
	n_data = length(x)
	
	KSs = Array{Float64,1}()
	for i in 1:n_sample
		@info(" $i/$n_sample")
		z = pl(rand(n_data),s0,C0,mi)
#		z = Array{Float64,1}()
#		for j in 1:n_data
#			push!(z,pl(rand(),s0,C0,mi))
#			push!(z,pl_trunc(rand(),s0,mi,ma))
#			push!(z,pl_trunc(rand(),s0,mi,n_data))
#		end
		s = ml_clauset(z)
		C = 1/sum((mi:1e7).^(-s))
		temp = KS_pl(z,s,C,mi)
		KS = temp[1]
#		KS = KS_pl(z,s0,C0,mi)
		push!(KSs,KS)
	end
	
	p = sum(KSs .> KS0)/n_sample
	
	return p
end



function old_pl(y::Float64,s::Float64,C::Float64,mi::Float64=1.)
	x = 0.
	n = mi - 1
	while x < y
		n += 1
		x += C*n^(-s)
	end
	
	return n
end

function pl_trunc(y::Float64,s::Float64,mi::Float64,ma::Float64)
	x = 0.
	n = mi - 1
	C = 1/sum((mi:ma).^(-s))
	while x < y
		n += 1
		x += C*n^(-s)
	end
	
	return n
end

function pl_trunc(y::Float64,s::Float64,mi::Float64,n_data::Int64)
	x = 0
	n = mi - 1
	C0 = 1/sum((1:1e5).^(-s))
	ma = (1/(C0*n_data))^(-1/s)
	C1 = 1/sum((mi:ma).^(-s))
	while x < y
		n += 1
		x += C1*n^(-s)
	end
	
	return n
end

function my_KS(x::Array{Float64,1},s::Float64,C::Float64,bins::Array{Float64,1},midbins::Array{Float64,1})
	mi = minimum(x)
	ma = maximum(x)
	n = length(mi:ma)
	
	th = C*(mi:ma).^(-s)
	bined_th = Array{Float64,1}()
	for i in 1:length(midbins)
		push!(bined_th,sum((th .> bins[i]).*(th .<= bins[i+1]).*th))
	end
	cdf = [bined_th[1],]
	for i in 2:length(bined_th)
		push!(cdf,cdf[end]+bined_th[i])
	end
	
	ecdf = Array{Float64,1}()
	for i in midbins
		push!(ecdf,sum(x .<= i)/length(x))
	end
	
	KS = maximum(abs(cdf-ecdf))
	
	return KS
end


function KS_pl_old(x::Array{Float64,1},s::Float64,C::Float64)
	mi = minimum(x)
	ma = maximum(x)
	n = length(mi:ma)
	me = ceil(Int64,.2*n)
	tail = ceil(Int64,.8*n)

	cdf = [C*mi^(-s),]
	for i in (mi+1):ma
		push!(cdf,cdf[end]+C*i^(-s))
	end
		
	ecdf = Array{Float64,1}()
	for i in mi:ma
		push!(ecdf,sum(x .<= i)/length(x))
	end
	
	KS = maximum(abs(cdf[me:end]-ecdf[me:end]))
	tail_diff_max = maximum(cdf[tail:end]-ecdf[tail:end])
	tail_diff_min = minimum(cdf[tail:end]-ecdf[tail:end])
	tail_diff_2 = norm(cdf[tail:end]-ecdf[tail:end],2)

	return KS,tail_diff_max,tail_diff_min,tail_diff_2
end



function gof_pl(j::String,x::Array{Float64,1},s0::Float64,C0::Float64,mi::Float64,n_sample::Int=2500)
	KS0 = KS_pl(x,s0,C0)
	n_data = length(x)

	KSs = Array{Float64,1}()
	for i in 1:n_sample
#		if (i%100 == 0) || i == n_sample
			@info(j*", GoF power law: $i/$n_sample")
#		end
		z = rand_pl(rand(n_data),s0,C0,mi)
		s = mle_pl(z)
		C = 1/zeta(s,mi)
		KS = KS_pl(z,s,C)
		push!(KSs,KS)
	end

	p = sum(KSs .> KS0)/n_sample
	
	writedlm("./analysis/"*j*"_KS0_pl.csv",KS0,',')
	writedlm("./analysis/"*j*"_KSs_pl.csv",KSs,',')
	
	return p
end

function gof_plc(j::String,x::Array{Float64,1},a0::Float64,l0::Float64,C0::Float64,mi::Float64,n_sample::Int=2500)
	KS0 = KS_plc(x,a0,l0,C0)
	n_data = length(x)
	
	KSs = Array{Float64,1}()
	for i in 1:n_sample
		if (i%100 == 0) || i == n_sample
			@info(j*", GoF power law with cutoff: $i/$n_sample")
		end
		z = rand_plc(rand(n_data),a0,l0,C0,mi)
		a,l = mle_plc(z)
		C = 1/(real(polylog(a,Complex(exp(-l)))) - sum((1:mi-1).^(-a).*exp.(-l*(1:mi-1))))
		KS = KS_plc(z,a,l,C)
		push!(KSs,KS)
	end

	p = sum(KSs .> KS0)/n_sample

	writedlm("./analysis/"*j*"_KS0_plc.csv",KS0,',')
	writedlm("./analysis/"*j*"_KSs_plc.csv",KSs,',')
	
	return p
end

function gof_yule(j::String,x::Array{Float64,1},a0::Float64,C0::Float64,mi::Float64,n_sample::Int=2500)
	KS0 = KS_yule(x,a0,C0)
	n_data = length(x)
	
	KSs = Array{Float64,1}()
	for i in 1:n_sample
		if (i%100 == 0) || i == n_sample
			@info(j*", GoF Yule law: $i/$n_sample")
		end
		z = rand_yule(rand(n_data),a0,C0,mi)
		a = mle_yule(z,mi)
		C = 1/(1-(a-1)*sum(beta.(1:(mi-1),a)))
		KS = KS_yule(z,a,C)
		push!(KSs,KS)
	end
	
	p = sum(KSs .> KS0)/n_sample
	
	writedlm("./analysis/"*j*"_KS0_yule.csv",KS0,',')
	writedlm("./analysis/"*j*"_KSs_yule.csv",KSs,',')
	
	return p
end

function gof_poisson(j::String,x::Array{Float64,1},mu0::Float64,C0::Float64,mi::Float64,n_sample::Int=2500) 
	KS0 = KS_poisson(x,mu0,C0)
	n_data = length(x)
	
	KSs = Array{Float64,1}()
	for i in 1:n_sample
		if (i%100 == 0) || i == n_sample
			@info(j*", GoF Poisson distri: $i/$n_sample")
		end
		z = rand_poisson(rand(n_data),mu0,C0,mi)
		mu = mle_poisson(z,mi)
		C = exp(mu) - sum(mu.^(0:mi-1)./(factorial.(0:mi-1)))
		KS = KS_poisson(z,mu,C)
		push!(KSs,KS)
	end
	
	p = sum(KSs .> KS0)/n_sample
	
	writedlm("./analysis/"*j*"_KS0_poisson.csv",KS0,',')
	writedlm("./analysis/"*j*"_KSs_poisson.csv",KSs,',')
	
	return p
end

function KS_pl(x::Array{Float64,1},s::Float64,C::Float64)
	mi = minimum(x)
	ma = maximum(x)
	n = length(x)

	cdf = [C*mi^(-s),]
	ecdf = [sum(x .<= mi)/n,]
	for i in (mi+1):ma
		push!(cdf,cdf[end]+C*i^(-s))
		push!(ecdf,sum(x .<= i)/n)
	end
	
#	KS = maximum(abs.(cdf[Array{Int,1}(mi:ma)] - ecdf[Array{Int,1}(mi:ma)]))
	KS = maximum(abs.(cdf - ecdf))
	
	return KS
end

function KS_plc(x::Array{Float64,1},a::Float64,l::Float64,C::Float64)
	mi = minimum(x)
	ma = maximum(x)
	n = length(x)

	cdf = [C*mi^(-a)*exp(-l*mi),]
	ecdf = [sum(x .<= mi)/n,]
	for i in (mi+1):ma
		push!(cdf,cdf[end]+C*i^(-a)*exp(-l*i))
		push!(ecdf,sum(x .<= i)/n)
	end

	KS = maximum(abs.(cdf - ecdf))

	return KS
end

function KS_yule(x::Array{Float64,1},a::Float64,C::Float64)
	mi = minimum(x)
	ma = maximum(x)
	n = length(x)
	
	cdf = [C*(a-1)*beta(mi,a),]
	ecdf = [sum(x .<= mi)/n,]
	for i in (mi+1):ma
		push!(cdf,cdf[end]+C*(a-1)*beta(i,a))
		push!(ecdf,sum(x .<= i)/n)
	end
	
	KS = maximum(abs.(cdf - ecdf))
	
	return KS
end

function KS_poisson(x::Array{Float64,1},mu::Float64,C::Float64)
	mi = minimum(x)
	ma = maximum(x)
	n = length(x)
	
	last_term = C*mu^mi/factorial(mi)
	cdf = [last_term,]
	ecdf = [sum(x .<= mi)/n,]
	for i in (mi+1):ma
		last_term *= mu/i
		push!(cdf,cdf[end] + last_term)
		push!(ecdf,sum(x .<= i)/n)
	end
	
	KS = maximum(abs.(cdf - ecdf))
	
	return KS
end

function rand_pl(y::Array{Float64,1},s::Float64,C::Float64,mi::Float64=1.)
	x = 0.
	n = mi - 1
	todo = trues(length(y))
	ns = n*ones(length(y))
	while x < maximum(y) #&& x < 1 - 1/(2500*100)
		ns += todo .* ones(length(y))
		n += 1
		x += C*n^(-s)
		todo = [y[i] > x for i in 1:length(y)]
	end

	return ns
end

function rand_plc(y::Array{Float64,1},a::Float64,l::Float64,C::Float64,mi::Float64=1.)
	x = 0.
	n = mi - 1
	todo = trues(length(y))
	ns = n*ones(length(y))
	while x < maximum(y) #&& x < 1 - 1/(2500*100)
		ns += todo .* ones(length(y)) 	
		n += 1
		x += C*n^(-a)*exp(-l*n)
		todo = [y[i] > x for i in 1:length(y)]
	end

	return ns
end

function rand_yule(y::Array{Float64,1},a::Float64,C::Float64,mi::Float64=1.)
	x = 0.
	n = mi - 1
	todo = trues(length(y))
	ns = n*ones(length(y))
	while x < maximum(y) #&& x < 1 - 1/(2500*100)
		ns += todo .* ones(length(y))
		n += 1
		x += C*(a-1)*beta(n,a)
		todo = [y[i] > x for i in 1:length(y)]
	end
	
	return ns
end

function rand_poisson(y::Array{Float64,1},mu::Float64,C::Float64,mi::Float64=1.)
	x = 0.
	n = mi - 1
	todo = trues(length(y))
	ns = n*ones(length(y))
	last_term = mu^n/factorial(n)
	while x < maximum(y) && 1 - 1/(2500*100)
		ns += todo .* ones(length(y))
		n += 1
		last_term *= mu/n
		x += last_term
		todo = [y[i] > x for i in 1:length(y)]
	end
	
	return ns
end

