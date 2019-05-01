using PyPlot, DelimitedFiles, SpecialFunctions, Logging, Dates

include("journals.jl")
include("my_histo.jl")
include("mle.jl")
include("gof.jl")

number_sample = 500

zipf_plot = false
tail_plot = false
plots = false
pl = true
	s = 0.
pl_co = true
	a = 0.
	l = 0.
expo = false # useless, tail is too weak
yule = true
	al = 0.
poisson = false # useless, tail is too weak
stretch_expo = false # not done, tail is too weak
lognormal = false # not done, is exactly a parabola when plotted in loglog scales

#js = ["energy",]
js = [used_journals[parse(Int,ARGS[1])],]
#js = [used_journals[jn],]

ps = Array{Float64,1}()
p_gof = Dict{String,Array{Float64,1}}

for j in js
	global s,a,l,al,p_gof,ps
# Plot histogram
	@info("$(now()) -- Collecting data: $j")
	num = Array{Float64,1}(vec(readdlm("./data/"*j*".txt",'\t')[2:end-2,2]))
	mi = minimum(num)
	ma = maximum(num)
	nbins = Int(ma)
	distributions = Array{String,1}()
	max_k = 0.
	
	if plots
		for i in mi:ma
			figure(j)
			PyPlot.plot([i,i],[1e-8,sum(num.==i)/length(num)],"-b",linewidth=2,color=journals_colors[j][1])
		end
	end
		
# Power law
	if pl
		@info("$(now()) -- $j: Power law...")	
## MLE
		s = mle_pl(num)
		C = 1/zeta(s,mi)
## Goodness-of-fit
		p_pl = gof_pl(j,num,s,C,mi,number_sample)
		push!(ps,p_pl)
		push!(distributions,"Power law")
		@info("$(now()) -- "*j*", power law: p-value = $p_pl")
## Plot
		if plots
			Hs = sum(1.0./((mi:ma).^s))
			z = C * ((mi:ma).^(-s))
			PyPlot.plot(mi:ma,z,"--k",label="Zipf's law fit, s = $(round(s; digits=3))",linewidth=2)
			max_k = ceil(Int64,(length(num)/Hs)^(1/s))
			PyPlot.plot([max_k,max_k],[.5/length(num),maximum(z)*length(num)*2],":k",linewidth=2)
		end
	end	
# Power law with cutoff
	if pl_co
		@info("$(now()) -- Power law with cutoff...")
## MLE
		a,l = mle_plc(num)
		C = 1/(real(polylog(a,Complex(exp(-l)))) - sum((1:mi-1).^(-a).*exp.(-l*(1:mi-1))))
## Goodness-of-fit
		p_plc = gof_plc(j,num,a,l,C,mi,number_sample)
		push!(ps,p_plc)
		push!(distributions,"Power law with cutoff")
		@info("$(now()) -- "*j*", power law with cutoff: p-value = $p_plc")
## Plot
		if plots
			zz = C .* (mi:ma).^(-a) .* exp.(-l.*(mi:ma))
			PyPlot.plot(mi:ma,zz,".-k",label="PL with cutoff",linewidth=3)
		end
	end

# MLE of exponential distribution
	if expo
		@info("$(now()) -- Exponential distribution...")
		la = mle_exp(num,mi)
		C = (1-exp(-la))*exp(la*mi)
		zzz = C .* exp.(-la.*(mi:ma))
		PyPlot.plot(mi:ma,zzz,":k",label="Exp distribution",linewidth=3)
	end

# Yule distribution
	if yule
		@info("$(now()) -- Yule distribution...")
## MLE
		al = mle_yule(num,mi)
		C = 1/(1-(al-1)*sum(beta.(1:(mi-1),al)))
## Goodness-of-fit
		p_yule = gof_yule(j,num,al,C,mi,number_sample)
		push!(ps,p_yule)
		push!(distributions,"Yule-Simon distribution")
		@info("$(now()) -- "*j*", Yule-Simon distribution: p-value = $p_yule")
## Plot
		if plots
			zzzz = C*(al-1)*beta.(mi:ma,al)
			PyPlot.plot(mi:ma,zzzz,"--b",label="Yule distribution",linewidth=3)
		end
	end

# Poisson distribution
	if poisson
		@info("$(now()) -- Poisson distribution...")
## MLE
		mu = mle_poisson(num,mi)
		C = exp(mu) - sum(mu.^(0:mi-1)./(factorial.(0:mi-1)))
## Goodness-of-fit
		p_poisson = gof_poisson(j,num,mu,C,mi,number_sample)
		push!(ps,p_poisson)
		push!(distributions,"Poisson distribution")
		@info("$(now()) -- "*j*", Poisson distribution: p-value = $p_poisson")
## Plot
		if plots
			z5 = C .* mu.^(mi:ma)./(factorial.(mi:ma))
			PyPlot.plot(mi:ma,z5,"--g",label="Poisson distribution",linewidth=3)
		end
	end

	if plots
		title(journals_code[j]*", max. # articles predicted: $(ceil(Int,max_k))")
		xlabel("Number of articles")
		ylabel("Number of authors")
		axis([.9,max(ma*2,2*max_k),.5/length(num),2.])
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

j = js[1]

writedlm("./analysis/"*j*"_params_$(number_sample)_$iter.csv",[s,a,l,al],',')
writedlm("./analysis/"*j*"_p-gof_$(number_sample)_$iter.csv",ps,',')

