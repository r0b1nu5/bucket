using PyPlot, DelimitedFiles, SpecialFunctions, Dates

include("journals.jl")
include("my_histo.jl")
include("mle.jl")
include("gof.jl")
include("distributions.jl")

tups = Array{Tuple{String,Int64,Bool,Bool,Bool,Bool,Bool,Bool,Int64},1}()

for j in journals_short
	for id in 26:50
		push!(tups,(j,100,true,false,false,true,false,false,id))
	end
end

function journal_analysis_parallel(tup::Tuple{String,Int64,Bool,Bool,Bool,Bool,Bool,Bool,Int64})
	j = tup[1]
	number_sample = tup[2]
	save = tup[3]
	pl = tup[4]
	pl_co = tup[5]
	yule = tup[6]
	expo = tup[7]
	poisson = tup[8]
	id = tup[9]
	
################ COLLECTING DATA ################################
	@info("$(now()) -- Collecting data: $j")
	if j == "prl" || j == "prd"
		num = Int.(readdlm("data/"*j*"_reduced_parsed.csv",','))
	else
		num = Int.(readdlm("data/"*j*"_parsed.csv",','))
	end
	
	mi = num[1,1]
	ma = num[1,end]
	n_data = sum(num[2,:])
	distributions = Array{String,1}()
	max_k = 0.

	s = -1000.
	p_pl = -1000.
	a = -1000.
	l = -1000.
	p_plco = -1000.
	al = -1000.
	p_yule = -1000.
	b = -1000.
	p_exp = -1000.
	mu = -1000.
	p_poisson = -1000.
	
############### POWER LAW ###################################
	if pl
		@info("$(now()) -- $j: Power law...")	

# ============= mle =====================================
		s = new_mle_pl(num)
		C = C_pl(s,mi,ma)

# ============= goodness-of-fit =========================
		p_pl = new_gof_pl(j,num,s,C,mi,number_sample)
		@info("$(now()) -- "*j*", power law: p-value = $p_pl")
		@info "==========================================="
		
		if save
			writedlm("./analysis/"*j*"_pl_params_$(number_sample)_$id.csv",s,',')
			writedlm("./analysis/"*j*"_pl_p-gof_$(number_sample)_$id.csv",p_pl,',')
		end
		
	end	


################# POWER LAW WITH CUTOFF ########################
	if pl_co
		@info("$(now()) -- Power law with cutoff...")

# =========== mle ==========================================
		a,l = new_mle_plc(num)
		C = C_plc(a,l,mi,ma)

# =========== goodness-of-fit =============================
		p_plc = new_gof_plc(j,num,a,l,C,mi,number_sample)
		@info("$(now()) -- "*j*", power law with cutoff: p-value = $p_plc")
		@info "==========================================="
		
		if save
			writedlm("./analysis/"*j*"_plc_params_$(number_sample)_$id.csv",[a,l],',')
			writedlm("./analysis/"*j*"_plc_p-gof_$(number_sample)_$id.csv",p_plc,',')
		end
	end


################## YULE #####################################
	if yule
		@info("$(now()) -- Yule distribution...")
# =========== mle ==========================================
		al = new_mle_yule(num)
		C = C_ys(al,mi,ma)

# =========== goodness-of-fit =============================
		p_yule = new_gof_yule(j,num,al,C,mi,number_sample)
		@info("$(now()) -- "*j*", Yule-Simon distribution: p-value = $p_yule")
		@info "==========================================="

		if save
			writedlm("./analysis/"*j*"_yule_params_$(number_sample)_$id.csv",al,',')
			writedlm("./analysis/"*j*"_yule_p-gof_$(number_sample)_$id.csv",p_yule,',')
		end
	end


################# EXPONENTIAL ###################################
	if expo
		@info "$(now()) -- Exponential distribution..."
# =========== mle ==========================================
		b = new_mle_exp(num,mi)
		C = C_exp(b,mi,ma)
#		C = (1 - exp(-b))/exp(-b*mi)
# =========== goodness-of-fit =============================
		p_exp = new_gof_exp(j,num,b,C,mi,number_sample)
		
		@info "$(now()) -- "*j*", exponential distribution: p-value = $p_exp"
		@info "==========================================="
		
		if save
			writedlm("./analysis/"*j*"_expo_params_$(number_sample)_$id.csv",b,',')
			writedlm("./analysis/"*j*"_expo_p-gof_$(number_sample)_$id.csv",p_exp,',')
		end
	end


################ POISSON ########################################
	if poisson
		@info "$(now()) -- Poisson distribution..."
# =========== mle ==========================================
		mu = new_mle_poisson(num,mi)
		C = 1/(exp(mu) - sum((mu.^(0:mi-1))./(factorial.(0:mi-1))))
# =========== goodness-of-fit =============================
		p_poisson = new_gof_poisson(j,num,mu,C,mi,number_sample)
		
		@info " $(now()) -- "*j*", Poisson distribution: p-value = $p_poisson"
		@info "==========================================="

		if save
			writedlm("./analysis/"*j*"_poisson_params_$(number_sample)_$id.csv",mu,',')
			writedlm("./analysis/"*j*"_poisson_p-gof_$(number_sample)_$id.csv",p_poisson,',')
		end
	end

	
	return (pl = (s = s, p = p_pl), 
		plco = (a = a, l = l, p = p_plco),
		yule = (al = al, p = p_yule),
		expo = (b = b, p = p_exp),
		poisson = (mu = mu, p = p_poisson))
end



