using SparseArrays, DelimitedFiles, Statistics

include("line_dist.jl")

#ntw = "euroroad_red"
#ρ = .1
#ntw = "ba1"
#ρ = .01
#ntw = "pegase1354"
#ρ = .1
ntw = "ws1"
ρ = .1

dosimu = false
doplot = true

n_t = 2
n_i = 2

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

	τ1 = λ2/100
	τ2 = λn*10
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

	for i = 1:n_i
		for t = 1:n_t
			@info "=========== iter=$i/$(n_i), t=$t/$(n_t) ============"
	
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

	 #=
	writedlm("temp/"*ntw*"_ts.csv",τs,',')
	writedlm("temp/"*ntw*"_eff_t.csv",eff_θ,',')
	writedlm("temp/"*ntw*"_eff_p.csv",eff_ψ,',')
	writedlm("temp/"*ntw*"_confi_t.csv",confi_t,',')
	writedlm("temp/"*ntw*"_confj_t.csv",confj_t,',')
	writedlm("temp/"*ntw*"_confi_p.csv",confi_p,',')
	writedlm("temp/"*ntw*"_confj_p.csv",confj_p,',')
	writedlm("temp/"*ntw*"_erri_t.csv",erri_t,',')
	writedlm("temp/"*ntw*"_errj_t.csv",errj_t,',')
	writedlm("temp/"*ntw*"_erri_p.csv",erri_p,',')
	writedlm("temp/"*ntw*"_errj_p.csv",errj_p,',')
	# =#
end

if doplot
	# #=
	τs = vec(readdlm("temp/"*ntw*"_ts.csv",','))
	eff_θ = readdlm("temp/"*ntw*"_eff_t.csv",',')
	eff_ψ = readdlm("temp/"*ntw*"_eff_p.csv",',')
	n_i,n_t = size(eff_θ)
	# =#
	# #=
	figure(ntw)
	PyPlot.plot(τs,vec(sum(eff_θ .> .9,dims=1))./n_i,"o",color="C0")
	PyPlot.plot(τs,vec(sum(eff_θ,dims=1)./n_i),"--",color="C0")
	PyPlot.plot(τs,vec(sum(eff_ψ .> .9,dims=1))./n_i,"o",color="C1")
	PyPlot.plot(τs,vec(sum(eff_ψ,dims=1))./n_i,"--",color="C1")
	xlabel("τ")
	ylabel("success rate")
	# =#
end
		





