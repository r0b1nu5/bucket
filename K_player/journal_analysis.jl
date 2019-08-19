using PyPlot, DelimitedFiles, SpecialFunctions, Dates

include("journals.jl")
include("my_histo.jl")
include("mle.jl")
include("gof.jl")

tups = Array{Tuple{String,Int64,Bool,Bool,Bool,Bool,Int64},1}()
for j in journals_short
	for id in 1:25
		push!(tups,(j,100,true,true,true,true,id))
	end
end

function journal_analysis_parallel(tup::Tuple{String,Int64,Bool,Bool,Bool,Bool,Int64})
	j = tup[1]
	number_sample = tup[2]
	save = tup[3]
	pl = tup[4]
	pl_co = tup[5]
	yule = tup[6]
	id = tup[7]
	
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

	s = 0.
	a = 0.
	l = 0.
	al = 0.
	
############### POWER LAW ###################################
	if pl
		@info("$(now()) -- $j: Power law...")	

# ============= mle =====================================
		s = new_mle_pl(num)
		C = 1/zeta(s,mi)

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
		C = 1/(real(polylog(a,Complex(exp(-l)))) - sum((1:mi-1).^(-a).*exp.(-l*(1:mi-1))))

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
		C = 1/(1-(al-1)*sum(beta.(1:(mi-1),al)))

# =========== goodness-of-fit =============================
		p_yule = new_gof_yule(j,num,al,C,mi,number_sample)
		push!(ps,p_yule)
		push!(distributions,"Yule-Simon distribution")
		@info("$(now()) -- "*j*", Yule-Simon distribution: p-value = $p_yule")

	@info "==========================================="
	end

############### SAVING ANALYSIS ################################
	if save
		writedlm("./analysis/"*j*"_params_$(number_sample)_$id.csv",[s,a,l,al],',')
		writedlm("./analysis/"*j*"_p-gof_$(number_sample)_$id.csv",ps,',')
	end
	
	return s,a,l,al,ps
end



