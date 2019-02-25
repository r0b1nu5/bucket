using PyPlot 

include("ml_zipf.jl")

function LinRange(a::Float64,b::Float64,n::Int64)
	return linspace(a,b,n)
end

include("journals.jl")
include("my_histo.jl")

nbins = 300
logbins = true
linbins = false
zipf_plot = false
tail_plot = true

js = ["prl","prd","pre"]
#js = used_journals
#js = ["chaos","prl","lancet","pre","science","nature"]
#js = keys(journals_full)
#js = ["chaos","chaos_91-99","chaos_95-03","chaos_00-08","chaos_04-12","chaos_09-17","prl","prl_58-77","prl_68-87","prl_78-97","prl_88-07","prl_98-17","scientometrics","science_noano","nature_noano"]
#js = ["lancet_noano","lancet_letters_noano","bmj_noano","neng_j_med_noano","neng_j_med_letters_noano","jama_noano"]

a1 = 0;a2 = 0;b1 = 0;b2 = 0

for j in js
# Plot histogram
	num = collect(Float64,readdlm("./data/"*j*".txt",'\t')[2:end-2,2])
	mi = minimum(num)
	ma = maximum(num)
	nbins = Int(ma)
if logbins

# #=	
	for i in mi:ma
		figure(j)
		PyPlot.plot([i,i],[1e-8,sum(num.==i)/length(num)],"-b",linewidth=2,color=journals_colors[j][1])
	end
# =#
 
#=	
	bins1 = exp(collect(LinRange(log(mi-.5),log(ma+.5),nbins+1)))
	midbins1 = sqrt(bins1[2:end].*bins1[1:end-1])
	h1 = histo(num,bins1)
	for i in 1:length(h1)
		figure(j*" 1")
		PyPlot.plot([bins1[i];bins1[i+1];bins1[i+1];bins1[i];bins1[i]],[0;0;h1[i];h1[i];0].+1e-8,"-",linewidth=2,color=journals_colors[j][1])
	end
# =#

 #= Plot least-square fit
	figure(j*" 1")
	x = .5*(bins1[1:end-1]+bins1[2:end])
	lx = Array{Float64,1}()
	ly = Array{Float64,1}()
	for i in 1:nbins
		if h1[i] > 0
			push!(lx,log(x[i]))
			push!(ly,log(h1[i]))
		end
	end
	M = [lx'*lx ones(length(lx))'*lx;ones(length(lx))'*lx length(lx)]
	(a1,b1) = inv(M)*[lx'*ly;ones(length(ly))'*ly]
	PyPlot.plot([x;exp(-b1/a1)],[x;exp(-b1/a1)].^a1*exp(b1),"--c",label="y = $(round(exp(b1),2))*x^$(round(a1,2))",linewidth=2)
	
	maX = (exp(-b1)/2)^(1/a1)
	pred = exp(-b1/a1)
	
# =#
	
# Plot maximum likelihood fit for a Zipf's law
#=
	x = Array{Float64,1}()
	for nu in num
		push!(x,maximum(midbins1.*(bins1[1:end-1] .<= nu)))
	end
	s = ml_clauset(x)
=#
	s = ml_clauset(num)
	C = 1/sum((1:1e7).^(-s))
# Plot the power-law
	Hs = sum(1 ./((mi:ma).^s))
	z = C*((mi:ma).^(-s))
#	figure(j*" 0")
	PyPlot.plot(mi:ma,z,"--k",label="Zipf's law fit, s = $(round(s,3))",linewidth=2)
	max_k = ceil(Int64,(length(num)/Hs)^(1/s))
	PyPlot.plot([max_k,max_k],[.5/length(num),maximum(z)*length(num)*2],":k",linewidth=2)
	
	title(journals_code[j]*", max. # articles predicted: $(ceil(Int,max_k))")
	ylabel("Number of authors")
	xlabel("Number of articles")
	axis([.9,max(ma*2,2*max_k),.5/length(num),maximum(z)*2])
	loglog()
	legend()
	
end

if linbins
	bins2 = collect(LinRange(mi,ma,nbins+1))
	h2 = histo(num,bins2)/length(num)
	for i in 1:length(h1)
		figure(j*" 2")
		PyPlot.plot([bins2[i];bins2[i+1];bins2[i+1];bins2[i];bins2[i]],[0;0;h2[i];h2[i];0].+1e-8,"-r",linewidth=2)
	end

	figure(j*" 2")
	x = .5*(bins2[1:end-1]+bins2[2:end])
	lx = Array{Float64,1}()
	ly = Array{Float64,1}()
	for i in 1:nbins
		if h2[i] > 0
			push!(lx,log(x[i]))
			push!(ly,log(h2[i]))
		end
	end
	M = [lx'*lx ones(length(lx))'*lx;ones(length(lx))'*lx length(lx)]
	(a2,b2) = inv(M)*[lx'*ly;ones(length(ly))'*ly]
	PyPlot.plot([x;exp(-b2/a2)],[x;exp(-b2/a2)].^a2*exp(b2),"--m",label="y = $(round(exp(b2),2))*x^$(round(a2,2))")
	
	maX = (exp(-b2)/2)^(1/a2)
	pred = exp(-b2/a2)
	
# Plot maximum likelihood fit for a Zipf's law
	s = ml_zipf(num)
	Hs = sum(1 ./((mi:ma).^s))
	z = 1 ./((mi:ma).^s)/Hs
	PyPlot.plot(mi:ma,z,":sm",label="Zipf's law fit, s = $(round(s,3))",linewidth=2)
	max_k = ceil(Int64,length(num)^(1/s))
	PyPlot.plot([max_k,max_k],[.5/length(num),maximum(h1)*2],":k")

	title(journals_code[j]*", $(nbins) bins, max. # articles predicted: $(ceil(Int,pred))")
	ylabel("Number of authors")
	xlabel("Number of articles")
	axis([.4,max(ma*2,maX),.5/length(num),maximum(h2)*2])
	loglog()
	legend()
end

if zipf_plot
	vals = Array{Int64,1}()
	for i in mi:ma
		push!(vals,sum(num.==i))
	end
	
	figure(j*" 3")
	PyPlot.plot([mi,ma],[0,0],"--k")
#	PyPlot.plot(mi:ma,vals/sum(vals)-z,"ob")
	PyPlot.plot(mi:ma,abs(vals/sum(vals)-z),"oy",label="Zipf")
	if logbins
		PyPlot.plot(mi:ma,abs(vals/sum(vals)-(mi:ma).^a1*exp(b1)),"oc",label="Least-square fit to log-bins")
	end
	if linbins
		PyPlot.plot(mi:ma,abs(vals/sum(vals)-(mi:ma).^a2*exp(b2)),"om",label="Least-square fit to lin-bins")
	end
	semilogy()
	title(journals_code[j]*", goodness-of-fit")
	xlabel("Number of authors")
	ylabel("Error of Zipf's fit")
	legend()
end

if tail_plot
	s = ml_clauset(num)
	C = 1/sum((1:1e7).^(-s))
	print(C)
	cdf = [C*mi^(-s),]
	ecdf = [sum(num .== mi)/length(num),]
	for i in mi+1:ma
		push!(cdf,cdf[end] + C*i^(-s))
		push!(ecdf,sum(num .<= i)/length(num))
	end
	
	figure(j*" 666")
	PyPlot.plot(mi:ma,1+1e-8-ecdf,color=journals_colors[j][1],linewidth=2)
	PyPlot.plot(mi:ma,1+1e-8-cdf,"--k",linewidth=2)
	loglog()
	axis([mi,2*ma,min(1-maximum(ecdf[1:end-1]),1-maximum(cdf[1:end]))/2,1])
	title(journals_code[j]*": Tails comparison")
end
end








