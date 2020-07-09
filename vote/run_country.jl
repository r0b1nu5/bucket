using LinearAlgebra, Distributed, Dates

include("big_rand.jl")
include("scripts_mini.jl")
include("result.jl")

needed_paths = ["./temp_data/","./data/single/","./data/several/","./data/wta/"]
for path in needed_paths
	if !isdir(path)
		mkdir(path)
	end
end

 #=
NS = 9				# Number of states
N = rand(200:1000,NS)		# Population in each state
M = 2*round.(Int64,N/100) .+ 1	# Number of representative for each state
MU = .5*ones(NS)		# Mean opinions in each state
SI = .2*ones(NS)		# Opinions' standard deviation in each state
DE = .1*ones(NS)		# Opinions' bias in each state

w_ref = 0.1
emi = .3
eps = emi + .01
# =#

## Returns the natural opinions (X0), final opinions (X1), inverse dynamics matrix (LDis), and the list of state ids (state_id) for a country with NS states whose population is given by N and whose number of representatives is given by M. The distribution of natural opinions is bigaussian with mean \pm MU, standard deviation SI and horizontal shift DE. The communication distances considered are in epss.
function generate_country_shift(N::Array{Int64,1}, M::Array{Int64,1}, MU::Array{Float64,1}, DE::Array{Float64,1}, SI::Array{Float64,1}, eps::Float64)
	NS = length(N)
	global X0 = Array{Array{Float64,1},1}()
	for i in 1:NS
		maxd = 1000.
		while maxd > eps
			N1 = round(Int64,N[i]/2)
			N2 = N[i] - N1
			global x = big_rand(N1,-MU[i]+DE[i],SI[i],N2,MU[i]+DE[i],SI[i])
			maxd = maximum(x[2:end]-x[1:end-1])
		end	
		push!(X0,x)
	end
	
	X1 = Array{Array{Float64,1},1}()
	state_id = Array{Int64,1}()
	LDis = Array{Array{Float64,2},1}()
	for i in 1:NS
		A = Float64.((0 .< abs.(repeat(X0[i],1,N[i]) - repeat(X0[i]',N[i],1)) .< eps))
		d = vec(sum(A,dims=1))
		D = diagm(0 => d)
		L = D - A
		LpD = Symmetric(2*D - A)
	#	LpD = 2*D - A
		LDi = inv(LpD)*Diagonal(D)
	
		push!(X1,LDi*X0[i])
		push!(LDis,LDi)
		state_id = [state_id;i*ones(Int64,N[i])]
	end

	return X0, X1, LDis, state_id
end

## Returns the natural opinions (X0), final opinions (X1), inverse dynamics matrix (LDis), and the list of state ids (state_id) for a country with NS states whose population is given by N and whose number of representatives is given by M. The distribution of natural opinions is bigaussian with mean \pm MU, standard deviation SI. The population in the left peak is bias*N and in the right peak (1-bias)*N. The communication distances considered are in epss.
function generate_country_bias(N::Array{Int64,1}, M::Array{Int64,1}, MU::Array{Float64,1}, bias::Array{Float64,1}, SI::Array{Float64,1}, eps::Float64)
	NS = length(N)
	global X0 = Array{Array{Float64,1},1}()
	for i in 1:NS
		maxd = 1000.
		while maxd > eps
			N1 = round(Int64,bias[i]*N[i])
			N2 = N[i] - N1
			global x = big_rand(N1,-MU[i]+DE[i],SI[i],N2,MU[i]+DE[i],SI[i])
			maxd = maximum(x[2:end]-x[1:end-1])
		end	
		push!(X0,x)
	end
	
	X1 = Array{Array{Float64,1},1}()
	state_id = Array{Int64,1}()
	LDis = Array{Array{Float64,2},1}()
	for i in 1:NS
		A = Float64.((0 .< abs.(repeat(X0[i],1,N[i]) - repeat(X0[i]',N[i],1)) .< eps))
		d = vec(sum(A,dims=1))
		D = diagm(0 => d)
		L = D - A
		LpD = Symmetric(2*D - A)
	#	LpD = 2*D - A
		LDi = inv(LpD)*Diagonal(D)
	
		push!(X1,LDi*X0[i])
		push!(LDis,LDi)
		state_id = [state_id;i*ones(Int64,N[i])]
	end

	return X0, X1, LDis, state_id
end
function generate_country_par(abcdefghw::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Float64})
	amin,a,b,c,d,e,f,g,h,w = abcdefghw

	eps = readdlm("temp_data/eps_$amin.$b.csv",',')[1]
	N = Int.(vec(readdlm("temp_data/N_$amin.$c.csv",',')))
	M = Int.(vec(readdlm("temp_data/M_$amin.$d.csv",',')))
	MU = vec(readdlm("temp_data/MU_$amin.$e.csv",','))
	if f > 0
		DE = vec(readdlm("temp_data/DE_$amin.$f.csv",','))
	elseif g > 0
		bias = vec(readdlm("temp_data/bias_$amin.$g.csv",','))
	end
	SI = vec(readdlm("temp_data/SI_$amin.$h.csv",','))

	NS = length(N)
	if f > 0
		X0,X1,LDis,state_id = generate_country_shift(N,M,MU,DE,SI,eps)
	elseif g > 0
		X0,X1,LDis,state_id = generate_country_bias(N,M,MU,bias,SI,eps)
	end
	
	for i in 1:NS
		writedlm("temp_data/X0_$a.$b.$c.$d.$e.$f.$g.$h.$i.csv",X0[i],',')
		writedlm("temp_data/X1_$a.$b.$c.$d.$e.$f.$g.$h.$i.csv",X1[i],',')
		writedlm("temp_data/LDis_$a.$b.$c.$d.$e.$f.$g.$h.$i.csv",LDis[i],',')
	end
	writedlm("temp_data/stid_$a.$b.$c.$d.$e.$f.$g.$h.csv",state_id,',')
end

# ======================================================================================
# Single representative	
# Compute the effort needed to change the outcome of an election where the whole country elects a single representative (e.g., 2nd round of French elections).
# ======================================================================================
# #=

function effort_single(X0::Array{Array{Float64,1},1}, X1::Array{Array{Float64,1},1}, LDis::Array{Array{Float64,2},1}, state_id::Array{Int64,1}, NS::Int64, N::Array{Int64,1}, M::Array{Int64,1}, w_ref::Float64)
	x0 = Array{Float64,1}()
	x1 = Array{Float64,1}()
	for i in 1:NS
		x0 = [x0;X0[i]]
		x1 = [x1;X1[i]]
	end
	oc0 = outcome(x1)[1]
	oc1 = copy(oc0)
	w0 = -sign(oc0)*w_ref
	inf_order = mini_sort(x0,(w0 < 0))
	c1 = 0
	
	while oc0*oc1 > 0.
		c2 = 0
		while oc0*oc1 > 0. && c2 < length(x1)
			c1 += 1
			c2 += 1
				to_inf = inf_order[c2]
			state = state_id[to_inf]
			index = setdiff((1:length(x1)).*(state_id .== state),[0.])
			x1[index] += w0*LDis[state][:,to_inf-sum(N[1:state-1])]
			oc1 = outcome(x1)[1]
		end
	end
	
	xi1 = c1*abs(w0)

	return xi1
end

function effort_single_par(abcdefghw::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Float64})
	amin,a,b,c,d,e,f,g,h,w_ref = abcdefghw
	
#	@info "$(now()) -- Running $a.$b.$c.$d.$e.$f.$g.$h, single"

	N = Int64.(vec(readdlm("temp_data/N_$amin.$c.csv",',')))
	M = Int64.(vec(readdlm("temp_data/M_$amin.$d.csv",',')))
	NS = length(N)
	X0 = Array{Array{Float64,1},1}()
	X1 = Array{Array{Float64,1},1}()
	LDis = Array{Array{Float64,2},1}()
	for i in 1:NS
		push!(X0,vec(readdlm("temp_data/X0_$a.$b.$c.$d.$e.$f.$g.$h.$i.csv",',')))
		push!(X1,vec(readdlm("temp_data/X1_$a.$b.$c.$d.$e.$f.$g.$h.$i.csv",',')))
		push!(LDis,readdlm("temp_data/LDis_$a.$b.$c.$d.$e.$f.$g.$h.$i.csv",','))
	end
	state_id = Int64.(vec(readdlm("temp_data/stid_$a.$b.$c.$d.$e.$f.$g.$h.csv",',')))

	xi1 = effort_single(X0,X1,LDis,state_id,NS,N,M,w_ref)
	writedlm("data/single/xi1_$a.$b.$c.$d.$e.$f.$g.$h.csv",xi1,',')

#	@info "$(now()) -- $a.$b.$c.$d.$e.$f.$g.$h single is done."
end

# =#

# ======================================================================================
# Several representatives
# Compute the effort needed to change the outcome of an election where each state elects some number of representatives (given in M) (e.g., Conseil National in CH).
# ======================================================================================
# #=

function effort_several(X0::Array{Array{Float64,1},1}, X1::Array{Array{Float64,1},1}, LDis::Array{Array{Float64,2},1}, state_id::Array{Int64,1}, NS::Int64, N::Array{Int64,1}, M::Array{Int64,1}, w_ref::Float64)
	X2 = copy(X1)
	Mp0 = Array{Int64,1}()
	margMp0 = Array{Float64,1}()
	Mn0 = Array{Int64,1}()
	margMn0 = Array{Float64,1}()
	for i in 1:NS
		xxx = result_margin(X1[i],M[i])
		push!(Mp0,xxx[1])
		push!(margMp0,xxx[2])
		push!(Mn0,xxx[3])
		push!(margMn0,xxx[4])
	end
	
	oc0 = sum(Mp0) - sum(Mn0)
	oc1 = copy(oc0)
	Mp1 = copy(Mp0)
	Mn1 = copy(Mn0)
	margMp1 = copy(margMp0)
	margMn1 = copy(margMn0)
	w0 = -sign(oc0)*w_ref
	c1 = 0
	
	while oc0*oc1 > 0
		if w0 > 0
			state_order = Int64.(sortslices([margMn1 1:NS],dims=1)[:,2])
		else
			state_order = Int64.(sortslices([margMp1 1:NS],dims=1)[:,2])
		end
		s = state_order[1]
		mp1 = Mp1[s]
		mp2 = copy(mp1)
		inf_order = mini_sort(X0[s],(w0 < 0))
		while mp2 == mp1
			c2 = 0
			while mp2 == mp1 && c2 < length(X2[s])
				c1 += 1
				c2 += 1
				to_inf = inf_order[c2]
				X2[s] += w0*LDis[s][:,to_inf]
				mp2,mn2 = result(X2[s],M[s])
			end
		end
		Mp1[s],margMp1[s],Mn1[s],margMn1[s] = result_margin(X2[s],M[s])
		oc1 = sum(Mp1) - sum(Mn1)
	end

	xi2 = c1*abs(w0)

	return xi2
end

function effort_several_par(abcdefghw::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Float64})
	amin,a,b,c,d,e,f,g,h,w_ref = abcdefghw
	
#	@info "$(now()) -- Running $a.$b.$c.$d.$e.$f.$g.$h, several"

	N = Int64.(vec(readdlm("temp_data/N_$amin.$c.csv",',')))
	M = Int64.(vec(readdlm("temp_data/M_$amin.$d.csv",',')))
	NS = length(N)
	X0 = Array{Array{Float64,1},1}()
	X1 = Array{Array{Float64,1},1}()
	LDis = Array{Array{Float64,2},1}()
	for i in 1:NS
		push!(X0,vec(readdlm("temp_data/X0_$a.$b.$c.$d.$e.$f.$g.$h.$i.csv",',')))
		push!(X1,vec(readdlm("temp_data/X1_$a.$b.$c.$d.$e.$f.$g.$h.$i.csv",',')))
		push!(LDis,readdlm("temp_data/LDis_$a.$b.$c.$d.$e.$f.$g.$h.$i.csv",','))
	end
	state_id = Int64.(vec(readdlm("temp_data/stid_$a.$b.$c.$d.$e.$f.$g.$h.csv",',')))

	xi2 = effort_several(X0,X1,LDis,state_id,NS,N,M,w_ref)
	writedlm("data/several/xi2_$a.$b.$c.$d.$e.$f.$g.$h.csv",xi2,',')
	
#	@info "$(now()) -- $a.$b.$c.$d.$e.$f.$g.$h several is done."
end
# =#


# ======================================================================================
# Winner takes all
# Compute the effort needed to change the outcome of an election where the "winner takes all" rule applies in each state (e.g., US presidential election). 
# ======================================================================================
# #=

function effort_wta(X0::Array{Array{Float64,1},1}, X1::Array{Array{Float64,1},1}, LDis::Array{Array{Float64,2},1}, state_id::Array{Int64,1}, NS::Int64, N::Array{Int64,1}, M::Array{Int64,1}, w_ref::Float64)
	X3 = copy(X1)
	Mp0 = Array{Int64,1}()
	margMp0 = Array{Float64,1}()
	Mn0 = Array{Int64,1}()
	margMn0 = Array{Float64,1}()
	for i in 1:NS
		xxx = result_margin(X1[i],1)
		push!(Mp0,M[i]*xxx[1])
		push!(margMp0,xxx[2])
		push!(Mn0,M[i]*xxx[3])
		push!(margMn0,xxx[4])
	end
	
	oc0 = sum(Mp0) - sum(Mn0)
	oc1 = copy(oc0)
	Mp1 = copy(Mp0)
	Mn1 = copy(Mn0)
	margMp1 = copy(margMp0)
	margMn1 = copy(margMn0)
	w0 = -sign(oc0)*w_ref
	c1 = 0
	
	while oc0*oc1 > 0
		if w0 > 0
			state_order = Int64.(sortslices([M./margMn1 1:NS],dims=1,rev=true)[:,2])
		else
			state_order = Int64.(sortslices([M./margMp1 1:NS],dims=1,rev=true)[:,2])
		end
		s = state_order[1]
		mp1 = Mp1[s]
		mp2 = copy(mp1)
		inf_order = mini_sort(X0[s],(w0 < 0))
		while mp2 == mp1
			c2 = 0
			while mp2 == mp1 && c2 < length(X3[s])
				c1 += 1
				c2 += 1
				to_inf = inf_order[c2]
				X3[s] += w0*LDis[s][:,to_inf]
				r = result(X3[s],1)
				mp2 = M[s]*r[1]
				mn2 = M[s]*r[2]
			end
		end
	
		R = result_margin(X3[s],1)
		Mp1[s] = M[s]*R[1]
		margMp1[s] = R[2]
		Mn1[s] = M[s]*R[3]
		margMn1[s] = M[s]*R[4]
		
		oc1 = sum(Mp1) - sum(Mn1)
	end
	
	xi3 = c1*abs(w0)

	return xi3
end

function effort_wta_par(abcdefghw::Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Float64})
	amin,a,b,c,d,e,f,g,h,w_ref = abcdefghw

#	@info "$(now()) -- Running $a.$b.$c.$d.$e.$f.$g.$h, wta"

	N = Int64.(vec(readdlm("temp_data/N_$amin.$c.csv",',')))
	M = Int64.(vec(readdlm("temp_data/M_$amin.$d.csv",',')))
	NS = length(N)
	X0 = Array{Array{Float64,1},1}()
	X1 = Array{Array{Float64,1},1}()
	LDis = Array{Array{Float64,2},1}()
	for i in 1:NS
		push!(X0,vec(readdlm("temp_data/X0_$a.$b.$c.$d.$e.$f.$g.$h.$i.csv",',')))
		push!(X1,vec(readdlm("temp_data/X1_$a.$b.$c.$d.$e.$f.$g.$h.$i.csv",',')))
		push!(LDis,readdlm("temp_data/LDis_$a.$b.$c.$d.$e.$f.$g.$h.$i.csv",','))
	end
	state_id = Int64.(vec(readdlm("temp_data/stid_$a.$b.$c.$d.$e.$f.$g.$h.csv",',')))

	xi3 = effort_wta(X0,X1,LDis,state_id,NS,N,M,w_ref)
	writedlm("data/wta/xi3_$a.$b.$c.$d.$e.$f.$g.$h.csv",xi3,',')
	
#	@info "$(now()) -- $a.$b.$c.$d.$e.$f.$g.$h wta is done."
end
# =#

# ===========================================================================================
# Run the simulation...
# ===========================================================================================
# #=

function run_country(ids::Tuple{Int64,Int64}, epss::Array{Float64,1}, Ns::Array{Array{Int64,1},1}, Ms::Array{Array{Int64,1},1}, MUs::Array{Array{Float64,1},1}, DEs::Array{Array{Float64,1},1}, biass::Array{Array{Float64,1},1}, SIs::Array{Array{Float64,1},1}, w_ref::Float64=.1, clean_disk::Bool=true)
	data_ids = Array{Tuple{Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Int64,Float64},1}()

	for b in 1:length(epss)
		writedlm("temp_data/eps_$(ids[1]).$b.csv",epss[b],',')
	end
	for c in 1:length(Ns)
		writedlm("temp_data/N_$(ids[1]).$c.csv",Ns[c],',')
	end
	for d in 1:length(Ms)
		writedlm("temp_data/M_$(ids[1]).$d.csv",Ms[d],',')
	end
	for e in 1:length(MUs)
		writedlm("temp_data/MU_$(ids[1]).$e.csv",MUs[e],',')
	end
	for f in 1:length(DEs)
		writedlm("temp_data/DE_$(ids[1]).$f.csv",DEs[f],',')
	end
	for g in 1:length(biass)
		writedlm("temp_data/bias_$(ids[1]).$g.csv",biass[g],',')
	end
	for h in 1:length(SIs)
		writedlm("temp_data/SI_$(ids[1]).$h.csv",SIs[h],',')
	end
	
	for a in ids[1]:ids[2]
		for b in 1:length(epss)
			for c in 1:length(Ns)
				for d in 1:length(Ms)
					for e in 1:length(MUs)
						for h in 1:length(SIs)
							for f in 1:length(DEs)
								push!(data_ids,(ids[1],a,b,c,d,e,f,0,h,w_ref))
							end
							for g in 1:length(biass)
								push!(data_ids,(ids[1],a,b,c,d,e,0,g,h,w_ref))
							end
						end
					end
				end
			end
		end
	end

	@info "$(now()) -- Generate countries"
	pmap(generate_country_par,data_ids)

	@info "$(now()) -- Single representative"
	pmap(effort_single_par,data_ids)

	@info "$(now()) -- Several representative"
	pmap(effort_several_par,data_ids)

	@info "$(now()) -- Winner takes all"
	pmap(effort_wta_par,data_ids)

# #=
	if clean_disk
		for a in ids[1]:ids[2]
			for b in 1:length(epss)
				for c in 1:length(Ns)
					NS = length(vec(readdlm("temp_data/N_$(ids[1]).$c.csv",',')))
					for d in 1:length(Ms)
						for e in 1:length(MUs)
							for h in 1:length(SIs)
								for f in 1:length(DEs)
									for i in 1:NS
										rm("temp_data/X0_$a.$b.$c.$d.$e.$f.0.$h.$i.csv")
										rm("temp_data/X1_$a.$b.$c.$d.$e.$f.0.$h.$i.csv")
										rm("temp_data/LDis_$a.$b.$c.$d.$e.$f.0.$h.$i.csv")
									end
									rm("temp_data/stid_$a.$b.$c.$d.$e.$f.0.$h.csv")
								end
								for g in 1:length(biass)
									for i in 1:NS
										rm("temp_data/X0_$a.$b.$c.$d.$e.0.$g.$h.$i.csv")
										rm("temp_data/X1_$a.$b.$c.$d.$e.0.$g.$h.$i.csv")
										rm("temp_data/LDis_$a.$b.$c.$d.$e.0.$g.$h.$i.csv")
									end
									rm("temp_data/stid_$a.$b.$c.$d.$e.0.$g.$h.csv")
								end
							end
						end
					end
				end
			end
		end
		
		for b in 1:length(epss)
			rm("temp_data/eps_$(ids[1]).$b.csv")
		end
		for c in 1:length(Ns)
			rm("temp_data/N_$(ids[1]).$c.csv")
		end
		for d in 1:length(Ms)
			rm("temp_data/M_$(ids[1]).$d.csv")
		end
		for e in 1:length(MUs)
			rm("temp_data/MU_$(ids[1]).$e.csv")
		end
		for f in 1:length(DEs)
			rm("temp_data/DE_$(ids[1]).$f.csv")
		end
		for g in 1:length(biass)
			rm("temp_data/bias_$(ids[1]).$g.csv")
		end
		for h in 1:length(SIs)
			rm("temp_data/SI_$(ids[1]).$h.csv")
		end
	end
# =#

	@info "$(now()) -- Finally done."
end

									




