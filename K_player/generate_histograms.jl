using PyPlot, DelimitedFiles, SpecialFunctions, Dates

include("journals.jl")
include("my_histo.jl")
include("mle.jl")
include("gof.jl")
include("distributions.jl")


number_sample = 2500
precision = 1e-4

############## JOURNALS TO USE ########################################
#js = [used_journals[parse(Int,ARGS[1])],]
#js = [used_journals[jn],]
#js = ["bmj", "food_chem_tox", "medical_ped_onc", "pnas", "pre", "chaos", "ieee_trans_autom_control", "nature", "energy", "lancet", "neng_j_med", "science", "scientometrics"]
js = journals_short
#js = ["prl", "prd"]
#js = ["prl_reduced", "prd_reduced"]

################ DO / DON'T LIST #######################################
plots = true
save = false
	iter = 100
pl = false
	s = 0.
pl_co = false
	a = 0.
	l = 0.
yule = false
	al = 0.
expo = false
	b = 0.
poisson = false
	m = 0.

# ============= useless ==============================================
stretch_expo = false # not done, tail is too weak
lognormal = false # not done, no closed form for the normalization constant. Furthermore, it is exactly a parabola when plotted in loglog scales, thus does not match the shape of the data very well...

zipf_plot = false
tail_plot = false
# ===================================================================


ps = Array{Float64,1}()
p_gof = Dict{String,Array{Float64,1}}

for j in js
	global s,a,l,al,p_gof,ps

################ COLLECTING DATA ################################
	
	@info("$(now()) -- Collecting data: $j")
	if j == "prl_reduced" || j == "prd_reduced"
		num = Int.(readdlm("data/"*j*"_parsed.csv",','))
		j = j[1:3]
	else
		num = Int.(readdlm("data/"*j*"_parsed.csv",','))
	end
	
	mi = num[1,1]
	ma = num[1,end]
	n_data = sum(num[2,:])
	distributions = Array{String,1}()
	max_k = 0.
	
############## PLOT HISTOGRAM ##################################
	if plots
		for i in 1:size(num)[2]
			figure(j)
			PyPlot.plot([num[1,i],num[1,i]],[1e-8,num[2,i]/n_data],linewidth=2,color=journals_colors[j][1])
		end
# ============ power law ===============================
		s = new_mle_pl(num)
		C = C_pl(s,mi,ma)
		z = z_pl(s,mi,ma)
		PyPlot.plot(mi:ma,z,"--k",label="pl: s = $(round(s; digits=3))",linewidth=1)
		
		Hs = sum(1.0./((mi:ma).^s))
		max_k = ceil(Int64,(sum(num[2,:])/Hs)^(1/s))
		PyPlot.plot([max_k,max_k],[.5/sum(num[2,:]),maximum(z)*sum(num[2,:])*2],"--k",linewidth=1)
		
# ============ power law with cutoff ===============================
		a,l = new_mle_plc(num)
		C = C_plc(a,l,mi,ma)
#		C = 1/(real(polylog(a,Complex(exp(-l)))) - sum((1:mi-1).^(-a).*exp.(-l*(1:mi-1))))
		zz = z_plc(a,l,mi,ma)
#		zz = C .* (mi:ma).^(-a) .* exp.(-l.*(mi:ma))
		PyPlot.plot(mi:ma,zz,"-.k",label="plc: a = $(round(a; digits=3)), l = $(round(l; digits=3))",linewidth=1)
		
# ============ yule-simon ===============================
		al = new_mle_yule(num)
		C = C_ys(al,mi,ma)
#		C = 1/(1-(al-1)*sum(beta.(1:(mi-1),al)))
		zzzz = z_ys(al,mi,ma)
#		zzzz = C*(al-1)*beta.(mi:ma,al)
		PyPlot.plot(mi:ma,zzzz,":k",label="yule: al = $(round(al; digits=3))",linewidth=1)
		
# ============ exponential ===============================
		b = new_mle_exp(num,mi)
		C = C_exp(b,mi,ma)
#		C = (1 - exp(-b))/exp(-b*mi)
		z = z_exp(b,mi,ma)
#		z = C * exp.(-b*(mi:ma))
		PyPlot.plot(mi:ma,z,"--k",label="exp: b = $(round(b; digits=3))",linewidth=1)
#=
# ============ poisson ===============================
		m = new_mle_poisson(num,mi)
		C = 1/(exp(m) - sum((m.^(0:mi-1))./(factorial.(0:mi-1))))
		z = C * (m.^(mi:ma))./(factorial.(Array{Float64,1}(mi:ma)))
		PyPlot.plot(mi:ma,z,"--r",label="Poisson distribution2",linewidth=3)

=#
		@info "==========================================="
		
		title(journals_code[j]*", max. # articles predicted: $(ceil(Int,max_k))")
		xlabel("Number of articles")
		ylabel("Number of authors")
		axis([.9*mi,2*ma,.5/sum(num[2,:]),2.])
		loglog()
		legend()
	end
		
############### POWER LAW ###################################
	if pl
		@info("$(now()) -- $j: Power law...")	

# ============= mle =====================================
		s = new_mle_pl(num)
		C = C_pl(s,mi,ma)
#		C = 1/zeta(s,mi)

# ============= goodness-of-fit =========================
		p_pl = new_gof_pl(j,num,s,C,mi,number_sample)
		push!(ps,p_pl)
		push!(distributions,"Power law")
		@info("$(now()) -- "*j*", power law: p-value = $p_pl")

	@info "==========================================="
		
	end	


################# POWER LAW WITH CUTOFF ########################
	if pl_co
		@info("$(now()) -- Power law with cutoff...")

# =========== mle ==========================================
		a,l = new_mle_plc(num)
		C = C_plc(a,l,mi,ma)
#		C = 1/(real(polylog(a,Complex(exp(-l)))) - sum((1:mi-1).^(-a).*exp.(-l*(1:mi-1))))

# =========== goodness-of-fit =============================
		p_plc = new_gof_plc(j,num,a,l,C,mi,number_sample)
		push!(ps,p_plc)
		push!(distributions,"Power law with cutoff")
		@info("$(now()) -- "*j*", power law with cutoff: p-value = $p_plc")

	@info "==========================================="
	end


################## YULE #####################################
	if yule
		@info("$(now()) -- Yule distribution...")
# =========== mle ==========================================
		al = new_mle_yule(num,mi)
		C = C_ys(al,mi,ma)
#		C = 1/(1-(al-1)*sum(beta.(1:(mi-1),al)))

# =========== goodness-of-fit =============================
		p_yule = new_gof_yule(j,num,al,C,mi,number_sample)
		push!(ps,p_yule)
		push!(distributions,"Yule-Simon distribution")
		@info("$(now()) -- "*j*", Yule-Simon distribution: p-value = $p_yule")

	@info "==========================================="
	end

############### SAVING ANALYSIS ################################
if save
	writedlm("./analysis/"*j*"_params_$(number_sample)_$iter.csv",[s,a,l,al],',')
	writedlm("./analysis/"*j*"_p-gof_$(number_sample)_$iter.csv",ps,',')
end







############### USELESS DISTRIBUTIONS ############################################
# MLE of exponential distribution
	if expo
		@info("$(now()) -- Exponential distribution...")
		la = mle_exp(num,mi)
		C = (1-exp(-la))*exp(la*mi)
		zzz = C .* exp.(-la.*(mi:ma))
		PyPlot.plot(mi:ma,zzz,":k",label="Exp distribution",linewidth=3)
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
			z5 = C .* mu.^(mi:ma)./(factorial.(Array{Float64,1}(mi:ma)))
			PyPlot.plot(mi:ma,z5,"--g",label="Poisson distribution",linewidth=3)
		end
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


