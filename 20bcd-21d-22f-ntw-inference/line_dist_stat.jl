using SparseArrays, DelimitedFiles, Statistics

include("line_dist.jl")
include("tools.jl")

#ntw = "euroroad_red"
#ntw = "euroroad_red00"
#ρ = .1
#ntw = "ba0" # Each new node connected to 1 other.
#ntw = "ba00"
#ρ = .001
#ntw = "ba1" # Each new node connected to 2 others.
#ρ = .01
#ntw = "ba5"
#ρ = .01
#ntw = "pegase1354"
#ntw = "pegase135400"
#ρ = .1
#ntw = "ws1"
#ρ = .1
#ntw = "ws2"
#ntw = "ws200"
#ρ = .1

@info "Which network?"
ntw = readline()

dosimu = false
doloadgroup = true
doload = false
doplot = true
plotrank = true

n_t = 20
n_i = 100
n_c = 10

if dosimu

	Lsp = readdlm("ntw_data/"*ntw*"_lap_mat_sp.csv",',')
#	L = Matrix(sparse(Int64.(Lsp[:,1]),Int64.(Lsp[:,2]),Lsp[:,3]))
	L = sparse(Int64.(Lsp[:,1]),Int64.(Lsp[:,2]),Lsp[:,3])
	Ld = pinv(Matrix(L))
	B,w,Bt = L2B(L)
	n,m = size(B)
	
	P = ρ*rand(n)
	P .-= mean(P)

#	τ1 = .0001
#	τ2 = .5
#	τs = LinRange(τ1,τ2,n_t)
#	e1 = -4.
#	e2 = 1.
#	τs = (10).^(LinRange(e1,e2,n_t))

	iter = 20000
	while iter > 19000
		global θs0,iter = kuramoto_sine(ntw,L,P,Ld*P,1,0.,.1,0.,false,20000,1e-5,.1)
		global θ0 = θs0[:,end]
	end

	J = spzeros(n,n)
	for k in 1:size(Lsp)[1]
		i = Int64(Lsp[k,1])
		j = Int64(Lsp[k,2])
		if i < j
			J[i,j] = -Lsp[k,3]*cos(θ0[i]-θ0[j])
			J[j,i] = -Lsp[k,3]*cos(θ0[i]-θ0[j])
		end
	end
	J = J - spdiagm(0 => vec(sum(J,dims=2)))
	λs = eigvals(Matrix(J))
	λ2 = -λs[end-1]
	λn = -λs[1]

	τ1 = 10/λ2
	τ2 = 1/(10*λn)
	τs = exp.(LinRange(log(τ1),log(τ2),n_t))
		
	for c in 1:n_c
		if m/n_c < 200
			ls = c:n_c:m
		else
			dc = floor(Int64,m/100)
			ls = c:dc:m
		end
	
		eff_θ = zeros(n_i,n_t)
		eff_ψ = zeros(n_i,n_t)
		confi_t = zeros(n_i,n_t)
		confj_t = zeros(n_i,n_t)
		confi_p = zeros(n_i,n_t)
		confj_p = zeros(n_i,n_t)
		erri_t = zeros(n_i,n_t)
		errj_t = zeros(n_i,n_t)
		erri_p = zeros(n_i,n_t)
		errj_p = zeros(n_i,n_t)
	
		for i = 1:n_i
			for t = 1:n_t
				@info "=========== counter = $c, iter=$i/$(n_i), t=$t/$(n_t) ============"
	
				l = ls[i]
				sl = findmax(B[:,l])[2]
				tl = findmin(B[:,l])[2]
		
				τ = τs[t]
		
				θs,it = kuramoto_sine(ntw,L,P,θ0,l,1.,1/τ,0.,true,1000,-1e-5,min(.01,τ/10))
				nθ = [norm(θs[i,:] .- θs[i,1],Inf) for i in 1:n]
				lθ = Int64.(sortslices([nθ 1:n],dims=1,rev=true)[:,2])
				rθ = Int64.(sortslices([lθ 1:n],dims=1)[:,2])
				iθ = lθ[1]
				jθ = lθ[2]
				kθ = lθ[3]
				uθ = union([sl,tl],[iθ,jθ])
				if length(uθ) < 2
					@info "Weird..."
					eff_θ[i,t] = NaN
				elseif length(uθ) == 2
						eff_θ[i,t] = 1.
				elseif length(uθ) == 3
					eff_θ[i,t] = .5
				end
				confi_t[i,t] = 1 - nθ[kθ]/nθ[iθ]
				confj_t[i,t] = 1 - nθ[kθ]/nθ[jθ]
#				erri_t[i,t] = 1 - maximum(nθ[[sl,tl]])/nθ[iθ]
#				errj_t[i,t] = 1 - minimum(nθ[[sl,tl]])/nθ[jθ]
				erri_t[i,t] = sort(rθ[[sl,tl]])[1]
				errj_t[i,t] = sort(rθ[[sl,tl]])[2]
		
				ψs = L*θs
				nψ = [norm(ψs[i,:] .- ψs[i,1],Inf) for i in 1:n]
				lψ = Int64.(sortslices([nψ 1:n],dims=1,rev=true)[:,2])
				rψ = Int64.(sortslices([lψ 1:n],dims=1)[:,2])
				iψ = lψ[1]
				jψ = lψ[2]
				kψ = lψ[3]
				uψ = union([sl,tl],[iψ,jψ])
				if length(uψ) < 2
					@info "Weird..."
					eff_ψ[i,t] = NaN
				elseif length(uψ) == 2
					eff_ψ[i,t] = 1.
				elseif length(uψ) == 3
					eff_ψ[i,t] = .5
				end
				confi_p[i,t] = 1 - nψ[kψ]/nψ[iψ]
				confj_p[i,t] = 1 - nψ[kψ]/nψ[jψ]
#				erri_p[i,t] = 1 - maximum(nψ[[sl,tl]])/nψ[iψ]
#				errj_p[i,t] = 1 - minimum(nψ[[sl,tl]])/nψ[jψ]
				erri_p[i,t] = sort(rψ[[sl,tl]])[1]
				errj_p[i,t] = sort(rψ[[sl,tl]])[2]
			end
		end	

		# #=
		writedlm("temp/"*ntw*"_ts.csv",τs,',')
		writedlm("temp/"*ntw*"_eff_t_c$c.csv",eff_θ,',')
		writedlm("temp/"*ntw*"_eff_p_c$c.csv",eff_ψ,',')
		writedlm("temp/"*ntw*"_confi_t_c$c.csv",confi_t,',')
		writedlm("temp/"*ntw*"_confj_t_c$c.csv",confj_t,',')
		writedlm("temp/"*ntw*"_confi_p_c$c.csv",confi_p,',')
		writedlm("temp/"*ntw*"_confj_p_c$c.csv",confj_p,',')
		writedlm("temp/"*ntw*"_erri_t_c$c.csv",erri_t,',')
		writedlm("temp/"*ntw*"_errj_t_c$c.csv",errj_t,',')
		writedlm("temp/"*ntw*"_erri_p_c$c.csv",erri_p,',')
		writedlm("temp/"*ntw*"_errj_p_c$c.csv",errj_p,',')
		# =#
	end
end

if doloadgroup
#	τs = vec(readdlm("temp/"*ntw*"_ts.csv",','))
	τs = vec(readdlm("temp/"*ntw*"_ts_c.csv",','))
	n_t = length(τs)
	global eff_θ = Matrix{Float64}(undef,0,n_t)
	eff_ψ = Matrix{Float64}(undef,0,n_t)
	confi_θ = Matrix{Float64}(undef,0,n_t)
	confj_θ = Matrix{Float64}(undef,0,n_t)
	confi_ψ = Matrix{Float64}(undef,0,n_t)
	confj_ψ = Matrix{Float64}(undef,0,n_t)
	erri_θ = Matrix{Float64}(undef,0,n_t)
	errj_θ = Matrix{Float64}(undef,0,n_t)
	erri_ψ = Matrix{Float64}(undef,0,n_t)
	errj_ψ = Matrix{Float64}(undef,0,n_t)
	
	for c in 1:n_c
		global eff_θ = [eff_θ;readdlm("temp/"*ntw*"_eff_t_c$c.csv",',')]
		global eff_ψ = [eff_ψ;readdlm("temp/"*ntw*"_eff_p_c$c.csv",',')]
		global confi_θ = [confi_θ;readdlm("temp/"*ntw*"_confi_t_c$c.csv",',')]
		global confj_θ = [confj_θ;readdlm("temp/"*ntw*"_confj_t_c$c.csv",',')]
		global confi_ψ = [confi_ψ;readdlm("temp/"*ntw*"_confi_p_c$c.csv",',')]
		global confj_ψ = [confj_ψ;readdlm("temp/"*ntw*"_confj_p_c$c.csv",',')]
		global erri_θ = [erri_θ;readdlm("temp/"*ntw*"_erri_t_c$c.csv",',')]
		global errj_θ = [errj_θ;readdlm("temp/"*ntw*"_errj_t_c$c.csv",',')]
		global erri_ψ = [erri_ψ;readdlm("temp/"*ntw*"_erri_p_c$c.csv",',')]
		global errj_ψ = [errj_ψ;readdlm("temp/"*ntw*"_errj_p_c$c.csv",',')]
	end
	qs_ciθ = get_quartiles(confi_θ)
	qs_cjθ = get_quartiles(confj_θ)
	qs_ciψ = get_quartiles(confi_ψ)
	qs_cjψ = get_quartiles(confj_ψ)
	qs_eiθ = get_quartiles(erri_θ)
	qs_ejθ = get_quartiles(errj_θ)
	qs_eiψ = get_quartiles(erri_ψ)
	qs_ejψ = get_quartiles(errj_ψ)
	
	eff_tot = (confj_θ .<= confj_ψ).*eff_ψ + (confj_θ .> confj_ψ).*eff_θ

	n_i,n_t = size(eff_θ)
end

if doload
	# #=
	τs = vec(readdlm("temp/"*ntw*"_ts.csv",','))
	eff_θ = readdlm("temp/"*ntw*"_eff_t.csv",',')
	eff_ψ = readdlm("temp/"*ntw*"_eff_p.csv",',')
	confi_θ = readdlm("temp/"*ntw*"_confi_t.csv",',')
	qs_ciθ = get_quartiles(confi_θ)
	confj_θ = readdlm("temp/"*ntw*"_confj_t.csv",',')
	qs_cjθ = get_quartiles(confj_θ)
	confi_ψ = readdlm("temp/"*ntw*"_confi_p.csv",',')
	qs_ciψ = get_quartiles(confi_ψ)
	confj_ψ = readdlm("temp/"*ntw*"_confj_p.csv",',')
	qs_cjψ = get_quartiles(confj_ψ)
	erri_θ = readdlm("temp/"*ntw*"_erri_t.csv",',')
	qs_eiθ = get_quartiles(erri_θ)
	errj_θ = readdlm("temp/"*ntw*"_errj_t.csv",',')
	qs_ejθ = get_quartiles(errj_θ)
	erri_ψ = readdlm("temp/"*ntw*"_erri_p.csv",',')
	qs_eiψ = get_quartiles(erri_ψ)
	errj_ψ = readdlm("temp/"*ntw*"_errj_p.csv",',')
	qs_ejψ = get_quartiles(errj_ψ)
	n_i,n_t = size(eff_θ)

	eff_tot = (confj_θ .<= confj_ψ).*eff_ψ + (confj_θ .> confj_ψ).*eff_θ
	# =#
end

if doplot
	ω2 = 10/τs[1]
	ωn = 1/(10*τs[end])
	# #=
	figure("method_eff_"*ntw)
	PyPlot.plot([ω2,ω2],[0.,1.],"--k")
	PyPlot.plot([ωn,ωn],[0.,1.],"--k")
	PyPlot.plot(1 ./τs,vec(sum(eff_θ .> .9,dims=1))./n_i,"o",color="C0")
	PyPlot.plot(1 ./τs,vec(sum(eff_θ .> .2,dims=1))./n_i,"-",color="C0")
	PyPlot.plot(1 ./τs,vec(sum(eff_ψ .> .9,dims=1))./n_i,"o",color="C1")
	PyPlot.plot(1 ./τs,vec(sum(eff_ψ .> .2,dims=1))./n_i,"-",color="C1")
	PyPlot.plot(1 ./τs,vec(sum(eff_tot .> .9,dims=1))./n_i,"--k")
#	PyPlot.plot(1 ./τs,vec(sum(eff_tot .> .2,dims=1))./n_i,"--",color="C2")
	xlabel("ω")
	ylabel("success rate")

	figure("confidence_"*ntw)
	#plot_quart(qs_ciθ,τs,"C0")
	plot_quart(qs_cjθ,τs,"C0")
	#plot_quart(qs_ciψ,τs,"C1")
	plot_quart(qs_cjψ,τs,"C1")
	xlabel("ω")
	ylabel("confidence")

	figure("err_"*ntw)
	#plot_quart(qs_eiθ,τs,"C0")
	plot_quart(qs_ejθ,τs,"C0")
	#plot_quart(qs_eiψ,τs,"C1")
	plot_quart(qs_ejψ,τs,"C1")
	xlabel("ω")
	ylabel("error")
	# =#
	
	if plotrank
		aerrj = [mean(errj_ψ[:,t]) for t in 1:n_t]
		qerrj = [quantile(errj_ψ[:,t],.99) for t in 1:n_t]
		merrj = [maximum(errj_ψ[:,t]) for t in 1:n_t]

		figure("rank_"*ntw)
		PyPlot.plot([1/τs[1],1/τs[end]],[2.,2.],"--k")
		PyPlot.plot(1 ./τs,qerrj,"o",color="C1")
		PyPlot.plot(1 ./τs,merrj,"o",mfc="none",color="C1")
		xlabel("ω")
		ylabel("rank")
	end
end
		





