using SparseArrays, DelimitedFiles, Statistics

include("line_dist.jl")
include("tools.jl")

#ntw = "euroroad_red"
#ρ = .1
#ntw = "ba0" # Each new node connected to 1 other.
#ρ = .001
#ntw = "ba1" # Each new node connected to 2 others.
#ρ = .01
ntw = "pegase1354"
ρ = .1
#ntw = "ws1"
#ρ = .1
#ntw = "ws2"
#ρ = .1

dosimu = false
doloadgroup = false
doload = true
doplot = true

n_t = 2
n_i = 2
n_c = 1

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

	ls = rand(1:m,n_i)

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
	
	for c in 1:n_c
		for i = 1:n_i
			for t = 1:n_t
				@info "=========== counter = $c, iter=$i/$(n_i), t=$t/$(n_t) ============"
	
				l = ls[i]
				sl = findmax(B[:,l])[2]
				tl = findmin(B[:,l])[2]
		
				τ = τs[t]
		
				θs,it = kuramoto_sine(ntw,L,P,θ0,l,1.,1/τ,0.,true,1000,-1e-5,min(.01,τ/10))
				nθ = [norm(θs[i,:] .- θs[i,1],Inf) for i in 1:n]
				iθ = findmax(nθ)[2]
				jθ = findmax([nθ[1:iθ-1];0.;nθ[iθ+1:n]])[2]
				kθ = findmax([nθ[1:min(iθ,jθ)-1];0.;nθ[min(iθ,jθ)+1:max(iθ,jθ)-1];0.;nθ[max(iθ,jθ)+1:n]])[2]
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
				erri_t[i,t] = 1 - maximum(nθ[[sl,tl]])/nθ[iθ]
				errj_t[i,t] = 1 - minimum(nθ[[sl,tl]])/nθ[jθ]
		
				ψs = L*θs
				nψ = [norm(ψs[i,:] .- ψs[i,1],Inf) for i in 1:n]
				iψ = findmax(nψ)[2]
				jψ = findmax([nψ[1:iψ-1];0.;nψ[iψ+1:n]])[2]
				kψ = findmax([nψ[1:min(iψ,jψ)-1];0.;nψ[min(iψ,jψ)+1:max(iψ,jψ)-1];0.;nψ[max(iψ,jψ)+1:n]])[2]
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
				erri_p[i,t] = 1 - maximum(nψ[[sl,tl]])/nψ[iψ]
				errj_p[i,t] = 1 - minimum(nψ[[sl,tl]])/nψ[jψ]
			end
		end	

		# #=
		writedlm("temp/"*ntw*"_ts_c.csv",τs,',')
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
	τs = vec(readdlm("temp/"*ntw*"_ts.csv",','))
	n_t = length(τs)
	eff_θ = Matrix{Float64}(undef,0,n_t)
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
		eff_θ = [eff_θ;readdlm("temp/"*ntw*"_eff_t_c$c.csv",',')]
		eff_ψ = [eff_ψ;readdlm("temp/"*ntw*"_eff_p_c$c.csv",',')]
		confi_θ = [confi_θ;readdlm("temp/"*ntw*"_confi_t_c$c.csv",',')]
		confj_θ = [confj_θ;readdlm("temp/"*ntw*"_confj_t_c$c.csv",',')]
		confi_ψ = [confi_ψ;readdlm("temp/"*ntw*"_confi_p_c$c.csv",',')]
		confj_ψ = [confj_ψ;readdlm("temp/"*ntw*"_confj_p_c$c.csv",',')]
		erri_θ = [erri_θ;readdlm("temp/"*ntw*"_erri_t_c$c.csv",',')]
		errj_θ = [errj_θ;readdlm("temp/"*ntw*"_errj_t_c$c.csv",',')]
		erri_ψ = [erri_ψ;readdlm("temp/"*ntw*"_erri_p_c$c.csv",',')]
		errj_ψ = [errj_ψ;readdlm("temp/"*ntw*"_errj_p_c$c.csv",',')]
	end
	qs_ciθ = get_quartiles(confi_θ)
	qs_cjθ = get_quartiles(confj_θ)
	qs_ciψ = get_quartiles(confi_ψ)
	qs_cjψ = get_quartiles(confj_ψ)
	qs_eiθ = get_quartiles(erri_θ)
	qs_ejθ = get_quartiles(errj_θ)
	qs_eiψ = get_quartiles(erri_ψ)
	qs_ejψ = get_quartiles(errj_ψ)

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
	μ2 = τs[1]/10
	μn = 10*τs[end]
	# #=
	figure("method_eff_"*ntw)
	PyPlot.plot([μ2,μ2],[0.,1.],"--k")
	PyPlot.plot([μn,μn],[0.,1.],"--k")
	PyPlot.plot(τs,vec(sum(eff_θ .> .9,dims=1))./n_i,"o",color="C0")
	PyPlot.plot(τs,vec(sum(eff_θ .> .2,dims=1))./n_i,"-",color="C0")
	PyPlot.plot(τs,vec(sum(eff_ψ .> .9,dims=1))./n_i,"o",color="C1")
	PyPlot.plot(τs,vec(sum(eff_ψ .> .2,dims=1))./n_i,"-",color="C1")
	PyPlot.plot(τs,vec(sum(eff_tot .> .9,dims=1))./n_i,"--k")
#	PyPlot.plot(τs,vec(sum(eff_tot .> .2,dims=1))./n_i,"--",color="C2")
	xlabel("τ")
	ylabel("success rate")

	figure("confidence_"*ntw)
	#plot_quart(qs_ciθ,τs,"C0")
	plot_quart(qs_cjθ,τs,"C0")
	#plot_quart(qs_ciψ,τs,"C1")
	plot_quart(qs_cjψ,τs,"C1")
	xlabel("τ")
	ylabel("confidence")

	figure("err_"*ntw)
	#plot_quart(qs_eiθ,τs,"C0")
	plot_quart(qs_ejθ,τs,"C0")
	#plot_quart(qs_eiψ,τs,"C1")
	plot_quart(qs_ejψ,τs,"C1")
	xlabel("τ")
	ylabel("error")
	# =#
end
		





