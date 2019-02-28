
using PyPlot, DelimitedFiles, SpecialFunctions, Logging

include("ml_zipf_1.0.jl")
include("journals.jl")
include("my_histo.jl")
include("mle.jl")

logbins = true
zipf_plot = false
tail_plot = true
pl = true
pl_co = false
expo = false # useless, tail is too weak
yule = false
poisson = false # useless, tail is too weak
stretch_expo = false # not done, tail is too weak
# TODO log-normal...

js = ["prl","prd"]
#js = journals_short

for j in js
# Plot histogram
	@info("Collecting data")
	num = Array{Float64,1}(vec(readdlm("./data/"*j*".txt",'\t')[2:end-2,2]))
	mi = minimum(num)
	ma = maximum(num)
	nbins = Int(ma)
if logbins
	@info("Entering logbins")
	for i in mi:ma
		figure(j)
		PyPlot.plot([i,i],[1e-8,sum(num.==i)/length(num)],"-b",linewidth=2,color=journals_colors[j][1])
	end

# MLE of power law
	@info("Power law...")	
	s = ml_clauset(num)
	C = 1/zeta(s,mi)
#	C = 1/sum((1:1e7).^(-s))
	
# Plot the power-law
if pl
	Hs = sum(1.0./((mi:ma).^s))
	z = C * ((mi:ma).^(-s))
	PyPlot.plot(mi:ma,z,"--k",label="Zipf's law fit, s = $(round(s; digits=3))",linewidth=2)
	max_k = ceil(Int64,(length(num)/Hs)^(1/s))
	PyPlot.plot([max_k,max_k],[.5/length(num),maximum(z)*length(num)*2],":k",linewidth=2)
end	
# MLE of power law with cutoff
if pl_co
	@info("Power law with cutoff...")
	a,l = mle_pl_cutoff(num)
	C = 1/real(polylog(a,Complex(exp(-l))))
	zz = C .* (mi:ma).^(-a) .* exp.(-l.*(mi:ma))
	PyPlot.plot(mi:ma,zz,".-k",label="PL with cutoff",linewidth=3)
end

# MLE of exponential distribution
if expo
	@info("Exponential distribution...")
	l = mle_exp(num,mi)
	C = (1-exp(-l))*exp(l*mi)
	zzz = C .* exp.(-l.*(mi:ma))
	PyPlot.plot(mi:ma,zzz,":k",label="Exp distribution",linewidth=3)
end

# MLE of Yule distribution
if yule
	@info("Yule distribution...")
	a = mle_yule(num,mi)
#	C = (a-1)*gamma(mi+a-1)/gamma(mi)
#	zzzz = C .* gamma.((mi:ma))./gamma.((mi:ma).+a)
	zzzz = (a-1)*beta.(mi:ma,a)
	PyPlot.plot(mi:ma,zzzz,"--b",label="Yule distribution",linewidth=3)
end

# MLE of Poisson distribution
if poisson
	@info("Poisson distribution...")
	mu = mle_poisson(num,mi)
	C = exp(mu) - sum(mu.^(0:mi-1)./(factorial.(0:mi-1)))
	z5 = C .* mu.^(mi:ma)./(factorial.(mi:ma))
	PyPlot.plot(mi:ma,z5,"--g",label="Poisson distribution",linewidth=3)
end

	title(journals_code[j]*", max. # articles predicted: $(ceil(Int,max_k))")
	xlabel("Number of articles")
	ylabel("Number of authors")
	axis([.9,max(ma*2,2*max_k),.5/length(num),maximum(z)*2])
	loglog()
	legend()
end

# if linbins

# if zipf_plot

if tail_plot
	s = ml_clauset(num)
	C = 1/zeta(s,mi)
	cdf = [C*mi^(-s),]
	ecdf = [sum(num .== mi)/length(num),]
	for i in mi+1:ma
		push!(cdf,cdf[end] + C*i^(-s))
		push!(ecdf,sum(num .<= i)/length(num))
	end
	
	figure(j*" 666")
	PyPlot.plot(mi:ma,1+1e-8 .- ecdf,color=journals_colors[j][1],linewidth=2)
	PyPlot.plot(mi:ma,1+1e-8 .- cdf,"--k",linewidth=2)
	loglog()
	axis([mi,2*ma,min(1-maximum(ecdf[1:end-1]),1-maximum(cdf[1:end]))/2,1])
	title(journals_code[j]*": Tails comparison")
end
end







