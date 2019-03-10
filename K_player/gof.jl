include("mle.jl")
using SpecialFunctions 

## ======================= Goodness-of-fit estimates ================================

function gof_pl(x::Array{Float64,1},s0::Float64,C0::Float64,mi::Float64,n_sample::Int=2500)
	KS0 = KS_pl(x,s0,C0)
	n_data = length(x)

	KSs = Array{Float64,1}()
	for i in 1:n_sample
#		if (i%100 == 0)
			@info("GoF power law: $i/$n_sample")
#		end
		z = rand_pl(rand(n_data),s0,C0,mi)
		s = mle_pl(z)
		C = 1/zeta(s,mi)
		KS = KS_pl(z,s,C)
		push!(KSs,KS)
	end

	p = sum(KSs .> KS0)/n_sample

	return p
end

function gof_plc(x::Array{Float64,1},a0::Float64,l0::Float64,C0::Float64,mi::Float64,n_sample::Int=2500)
	KS0 = KS_plc(x,a0,l0,C0)
	n_data = length(x)
	
	KSs = Array{Float64,1}()
	for i in 1:n_sample
		if (i%100 == 0) || i == n_sample
			@info("GoF power law with cutoff: $i/$n_sample")
		end
		z = rand_plc(rand(n_data),a0,l0,C0,mi)
		a,l = mle_plc(z)
		C = 1/real(polylog(a,Complex(exp(-l))))
		KS = KS_plc(z,a,l,C)
		push!(KSs,KS)
	end

	p = sum(KSs .> KS0)/n_sample

	return p
end


## ============================ Kolmogorov-Smirnov measures ==================================

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
	
	KS = maximum(abs.(cdf[Array{Int,1}(mi:ma)] - ecdf[Array{Int,1}(mi:ma)]))

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

	KS = maximum(abs.(cdf[Array{Int,1}(mi:ma)] - ecdf[Array{Int,1}(mi:ma)]))

	return KS
end


## ================================== Synthetic variables generation ================================

function rand_pl(y::Array{Float64,1},s::Float64,C::Float64,mi::Float64=1.)
	x = 0.
	n = mi - 1
	todo = trues(length(y))
	ns = n*ones(length(y))
	while x < maximum(y)
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
	while x < maximum(y)
		ns += todo .* ones(length(y))
		n += 1
		x += C*n^(-a)*exp(-l*n)
		todo = [y[i] > x for i in 1:length(y)]
	end

	return ns
end

## Old functions

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



