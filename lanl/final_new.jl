using DelimitedFiles, PyPlot, LinearAlgebra, JuMP, Ipopt, Distributed, FFTW


function run_new(Xs::Array{Float64,2},tau::Float64,kmax::Int64,plot::Bool=false,mu::Float64=1e-1,bp::Float64=1e-1)
	nn,NN = size(Xs)
	n = Int(nn/2)
	N = NN-1

	x = Xs[:,1:end-1]
	Dx = (Xs[(n+1):(2*n),2:end] - Xs[(n+1):(2*n),1:end-1])/tau
	xt = Array{Complex{Float64},2}(undef,nn,N)
	for i in 1:nn
		xt[i,:] = ifft(x[i,:])*sqrt(N)
	end
	Dxt = Array{Complex{Float64},2}(undef,n,N)
	for i in 1:n
		Dxt[i,:] = ifft(Dx[i,:])*sqrt(N)
	end

	Ls_l0 = zeros(n,kmax)
	L_l0 = 1000.
	A_l0 = zeros(n,nn)
	d_l0 = zeros(n)
	gamma_l0 = 0.
	l_l0 = 0
	k_l0 = 0
	for l in 1:n
		for k in 1:kmax
			@info "l0: l = $l, k = $k"
			Lt = Lmin_l0(x,Dx,xt,Dxt,l,k,mu,bp)
			Ls_l0[l,k] = Lt[1]
			if Lt[1] < L_l0
				L_l0,A_l0,d_l0,gamma_l0 = Lt
				l_l0 = l
				k_l0 = k
			end
		end
	end

	@info "Best l0: l = $(l_l0), ω = $(2*pi*(k_l0-1)/(tau*N)), L = $(L_l0)."
	
	Ls_l2 = zeros(kmax)
	L_l2 = 1000.
	A_l2 = zeros(n,nn)
	d_l2 = zeros(n)
	gamma_l2 = zeros(n)
	k_l2 = 0
	l_l2 = 0
	for k in 1:kmax
		@info "l2: k = $k"
		Lt = Lmin_l2(x,Dx,xt,Dxt,k,mu,bp)
		Ls_l2[k] = Lt[1]
		if Lt[1] < L_l2
			L_l2,A_l2,d_l2,gamma_l2 = Lt
			k_l2 = k
			l_l2 = findmax(gamma_l2)[2]
		end
	end

	@info "Best l2: ω = $(2*pi*(k_l2-1)/(tau*N)), L = $(L_l2)."

	if plot
		ws = Array(2*pi*(0:kmax-1)/(tau*N))
		
		plot_new_l0(Ls_l0, L_l0, l_l0, k_l0, ws, n)
		plot_new_l2(Ls_l2, L_l2, gamma_l2, l_l2, k_l2, ws)
	end

	return Ls_l0, (L_l0,A_l0,d_l0,gamma_l0,k_l0,l_l0), Ls_l2, (L_l2,A_l2,d_l2,gamma_l2,k_l2,l_l2)
end


function Lmin_l0(x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64},2}, l::Int64, k::Int64, mu::Float64=1e-1, bp::Float64=1e-1)
	nn,N = size(x)
	n = Int(nn/2)

	Sigma0 = (x*x')/N
	Sigma1 = (x*Dx')/N

	xtk = xt[:,k]
	Dxtlk = Dxt[l,k]
	
	Fk = real.(xtk*xtk')
	flk = real.(Dxtlk*xtk')
	glk = norm(Dxtlk)^2

	system_id = Model(with_optimizer(Ipopt.Optimizer, mu_init = mu, bound_push = bp))

	@variable(system_id, A1[i = 1:n, j = 1:n])
	for i in 1:n-1
		for j in i+1:n
			@constraint(system_id, A1[i,j] == A1[j,i])
		end
	end
	@variable(system_id, a2[i = 1:n])
	for i in 1:n
		@constraint(system_id, a2[i] >= 0.)
	end
	@variable(system_id, gamma >= 0.)

	@NLexpression(system_id, AtAS0[i = 1:n], sum(A1[j1,i]*A1[j1,j2]*Sigma0[j2,i] for j1 = 1:n for j2 = 1:n) + sum(A1[j,i]*a2[j]*Sigma0[n+j,i] for j = 1:n))
	@NLexpression(system_id, AtAS00[i = (n+1):(2*n)], sum(a2[i-n]*A1[i-n,j]*Sigma0[j,i] for j = 1:n) + a2[i-n]*a2[i-n]*Sigma0[i,i])
	@NLexpression(system_id, T1, sum(AtAS0[i] for i in 1:n) + sum(AtAS00[i] for i = (n+1):(2*n)))

	@NLexpression(system_id, AS1[i = 1:n], sum(A1[i,j]*Sigma1[j,i] for j = 1:n) + a2[i]*Sigma1[n+i,i])
	@NLexpression(system_id, T2, sum(AS1[i] for i in 1:n))

	@NLexpression(system_id, AAF[i = 1:n], A1[l,i]*(sum(A1[l,j]*Fk[j,i] for j = 1:n) + a2[l]*Fk[n+l,i]))
	@NLexpression(system_id, AAFF, a2[l]*(sum(A1[l,j]*Fk[j,n+l] for j = 1:n) + a2[l]*Fk[n+l,n+l]))
	@NLexpression(system_id, T3, sum(AAF[i] for i = 1:n) + AAFF)

	@NLexpression(system_id, T4, sum(flk[i]*A1[l,i] for i in 1:n) + flk[n+l]*a2[l])

	@NLexpression(system_id, g2, gamma^2)

	@NLobjective(system_id, Min, T1 + 2*T2 + .5*g2 - 2*gamma/sqrt(N)*sqrt(T3 + 2*T4 + glk))

	optimize!(system_id)

	return objective_value(system_id), value.(A1), value.(a2), value(gamma)
end


function Lmin_l2(x::Array{Float64,2}, Dx::Array{Float64,2}, xt::Array{Complex{Float64},2}, Dxt::Array{Complex{Float64},2}, k::Int64, mu::Float64=1e-1, bp::Float64=1e-1)
	nn,N = size(x)
	n = Int(nn/2)

	Sigma0 = (x*x')/N
	Sigma1 = (x*Dx')/N

	xtk = xt[:,k]
	
	Fk = real.(xtk*xtk')
	flk = [real.(Dxt[l,k]*xtk') for l = 1:n]
	glk = [norm(Dxt[l,k])^2 for l = 1:n]

	system_id = Model(with_optimizer(Ipopt.Optimizer, mu_init = mu, bound_push = bp))

	@variable(system_id, A1[i = 1:n, j = 1:n])
	for i in 1:n-1
		for j in i+1:n
			@constraint(system_id, A1[i,j] == A1[j,i])
		end
	end
	@variable(system_id, a2[i = 1:n])
	@variable(system_id, gamma[i = 1:n])
	for i in 1:n
		@constraint(system_id, a2[i] >= 0.)
		@constraint(system_id, gamma[i] >= 0.)
	end

	@NLexpression(system_id, AtAS0[i = 1:n], sum(A1[j1,i]*A1[j1,j2]*Sigma0[j2,i] for j1 = 1:n for j2 = 1:n) + sum(A1[j,i]*a2[j]*Sigma0[n+j,i] for j = 1:n))
	@NLexpression(system_id, AtAS00[i = (n+1):(2*n)], sum(a2[i-n]*A1[i-n,j]*Sigma0[j,i] for j = 1:n) + a2[i-n]*a2[i-n]*Sigma0[i,i])
	@NLexpression(system_id, T1, sum(AtAS0[i] for i in 1:n) + sum(AtAS00[i] for i = (n+1):(2*n)))

	@NLexpression(system_id, AS1[i = 1:n], sum(A1[i,j]*Sigma1[j,i] for j = 1:n) + a2[i]*Sigma1[n+i,i])
	@NLexpression(system_id, T2, sum(AS1[i] for i in 1:n))

	@NLexpression(system_id, g2, sum(gamma[i]^2 for i = 1:n))

	@NLexpression(system_id, AAF[l = 1:n, i = 1:n], A1[l,i]*(sum(A1[l,j]*Fk[j,i] for j = 1:n) + a2[l]*Fk[n+l,i]))
	@NLexpression(system_id, AAFF[l = 1:n], a2[l]*(sum(A1[l,j]*Fk[j,n+l] for j = 1:n) + a2[l]*Fk[n+l,n+l]))
	@NLexpression(system_id, T3[l = 1:n], sum(AAF[l,i] for i = 1:n) + AAFF[l])

	@NLexpression(system_id, T4[l = 1:n], sum(flk[l][i]*A1[l,i] for i in 1:n) + flk[l][n+l]*a2[l])

	@NLobjective(system_id, Min, T1 + 2*T2 + .5*g2 - 2/sqrt(N)*sum(gamma[l]*sqrt(T3[l] + 2*T4[l] + glk[l]) for l = 1:n))

	optimize!(system_id)

	return objective_value(system_id), value.(A1), value.(a2), value.(gamma)
end

function objective(Xs::Array{Float64,2}, tau::Float64, A1::Array{Float64,2}, d::Array{Float64,1}, gamma::Float64, l::Int64, k::Int64)
	x = Xs[:,1:end-1]
	nn,N = size(x)
	n = Int(nn/2)

	Dx = (Xs[(n+1):(2*n),2:end] - Xs[(n+1):(2*n),1:end-1])/tau
	xt = Array{Complex{Float64},2}(undef,nn,N)
	for i in 1:nn
		xt[i,:] = ifft(x[i,:])*sqrt(N)
	end
	Dxt = Array{Complex{Float64},2}(undef,n,N)
	for i in 1:n
		Dxt[i,:] = ifft(Dx[i,:])*sqrt(N)
	end
	
	Sigma0 = (x*x')/N
	Sigma1 = (x*Dx')/N

	xtk = xt[:,k]
	Dxtlk = Dxt[l,k]
	
	Fk = real.(xtk*xtk')
	flk = real.(Dxtlk*xtk')
	glk = norm(Dxtlk)^2

	A = [A1 diagm(0 => d)]

	obj = tr(transpose(A)*A*Sigma0) + 2*tr(A*Sigma1) + .5*gamma^2 - 2*gamma/sqrt(N)*sqrt(tr(A[l,:]*transpose(A[l,:])*Fk) + (2*flk*A[l,:])[1] + glk)

	return obj
end


function plot_new_l0(Ls_l0::Array{Float64,2}, L_l0::Float64, l_l0::Int64, k_l0::Int64, ws::Array{Float64,1}, n::Int64)
	figure(121)
	for l in 1:n
		PyPlot.plot(ws,Ls_l0[l,:],"o",label="l = $l")
	end
	xlabel("ω")
	ylabel("Log-likelihood")
	title("Best l0: l = $(l_l0), ω = $(round(ws[k_l0],sigdigits=5)), L = $(round(L_l0,sigdigits=5))")
end

function plot_new_l2(Ls_l2::Array{Float64,1}, L_l2::Float64, gamma_l2::Array{Float64,1}, l_l2::Int64, k_l2::Int64, ws::Array{Float64,1})
	figure(212)
	subplot(1,2,1)
	PyPlot.plot(ws,Ls_l2,"o")
	xlabel("ω")
	ylabel("Log-likelihood")
	title("Best l2: ω = $(round(ws[k_l2],sigdigits=5)), L = $(round(L_l2,sigdigits=5))")
	subplot(1,2,2)
	PyPlot.plot(1:n,gamma_l2,"o")
	xlabel("node index")
	ylabel("Estimated amplitude")
	title("Source: l = $(l_l2), γ_l = $(round(gamma_l2[l_l2],sigdigits=4))")
end

function plot_new_error(Lh::Array{Array{Float64,2},1}, L::Array{Float64,2}, dh::Array{Array{Float64,1},1}, d::Array{Float64,1}, gammah::Array{Array{Float64,1},1}, gamma::Array{Float64,1}, wh::Array{Float64,1}, w::Float64, Xss::Array{Array{Float64,2},1}, taus::Array{Float64,1}, ls::Array{Int64,1}, ks::Array{Int64,1})
	n = length(Lh)
	
	bs = Array(LinRange(-.4,.4,n+1))
	xs = (bs[2:end] + bs[1:end-1])/2
	wi = bs[2] - bs[1]
	a = Array(LinRange(1.,.3,n+1))
	
	figure(333)
	subplot(1,3,1)
	for i in 1:n
		PyPlot.bar(1 + xs[i], maximum(abs.(Lh[i] - L)./max.(abs.(L),1e-8).*(abs.(L) .> 1e-8)), wi, color="C0", alpha=a[i])
		PyPlot.bar(2 + xs[i], maximum(abs.(dh[i] - d)./max.(abs.(d),1e-8).*(abs.(d) .> 1e-8)), wi, color="C1", alpha=a[i])
		PyPlot.bar(3 + xs[i], maximum(abs.(gammah[i] - gamma)./max.(abs.(gamma),1e-8).*(abs.(gamma) .> 1e-8)), wi, color="C2", alpha=a[i])
		PyPlot.bar(4 + xs[i], abs(wh[i] - w)/abs(w), wi, color="C3", alpha=a[i])
	end

	xticks([1,2,3,4],["L","d","γ","ω"])
	ylabel("Max. relative error")

	subplot(1,3,2)
	for i in 1:n
		PyPlot.bar(1 + xs[i], norm(Lh[i] - L)/norm(L), wi, color="C0", alpha=a[i])
		PyPlot.bar(2 + xs[i], norm(dh[i] - d)/norm(d), wi, color="C1", alpha=a[i])
		PyPlot.bar(3 + xs[i], norm(gammah[i] - gamma)/norm(gamma), wi, color="C2", alpha=a[i])
		PyPlot.bar(4 + xs[i], abs(wh[i] - w)/abs(w), wi, color="C3", alpha=a[i])
	end

	xticks([1,2,3,4],["L","d","γ","ω"])
	ylabel("Relative error")

	bs = Array(LinRange(-.4,.4,n+2))
	xs = (bs[2:end] + bs[1:end-1])/2
	wi = bs[2] - bs[1]

	o0s = [objective(Xss[i],taus[i],L,d,maximum(gamma),ls[i],ks[i]) for i in 1:n]
	os = [objective(Xss[i],taus[i],Lh[i],dh[i],maximum(gammah[i]),ls[i],ks[i]) for i in 1:n]
	
	subplot(1,3,3)
	PyPlot.plot(1:n,os-o0s,"-o",color="C4")

	xticks(1:n)
	ylabel("Log-likelihood")
end
