using PyPlot, LinearAlgebra, DelimitedFiles

function coev(x0::Vector{Float64}, a0::Matrix{Float64}, save::Bool=true, δmax::Float64=1., δmin::Float64=0., γ::Float64=1., σ::Float64=1., tol::Float64=1e-6, maxiter::Int64=10000, h::Float64=.001)
	n = length(x0)
	Δ = δmax - δmin

	xs = copy(x0)
	as = copy(a0)

	iter = 0
	c = 0
	err = 1000.

	while iter < maxiter && err > tol
		iter += 1

		if iter%100 == 0
			@info "iter: $iter, err: $(round(err,digits=3))"
			c += 1
			writedlm("temp/xs_$c.csv",xs[:,1:end-1],',')
			xs = xs[:,end]
			#writedlm("temp/as_$c.csv",as[:,1:end-n],',')
			as = as[:,end-n+1:end]
		end

		x = xs[:,end]
		dx = x*ones(1,n) - ones(n)*x'
		a = as[:,end-n+1:end]

		kx1 = a*x - vec(sum(diagm(0 => x)*a,dims=2)) + log.((1 .- x)./(x .+ 1))./10
		ka1 = (Δ*exp.(-σ*dx.^2) .+ δmin - a)/γ
		kd1 = kx1*ones(1,n) - ones(n)*kx1'
		kx2 = (a+h/2*ka1)*(x+h/2*kx1) - vec(sum(diagm(0 => (x+h/2*kx1))*(a+h/2*ka1),dims=2)) + log.((1 .- (x+h/2*kx1))./((x+h/2*kx1) .+ 1))./10
		ka2 = (Δ*exp.(-σ*(dx+h/2*kd1).^2) .+ δmin - (a+h/2*ka1))/γ
		kd2 = kx2*ones(1,n) - ones(n)*kx2'
		kx3 = (a+h/2*ka2)*(x+h/2*kx2) - vec(sum(diagm(0 => (x+h/2*kx2))*(a+h/2*ka2),dims=2)) + log.((1 .- (x+h/2*kx2))./((x+h/2*kx2) .+ 1))./10
		ka3 = (Δ*exp.(-σ*(dx+h/2*kd2).^2) .+ δmin - (a+h/2*ka2))/γ
		kd3 = kx3*ones(1,n) - ones(n)*kx3'
		kx4 = (a+h*ka3)*(x+h*kx3) - vec(sum(diagm(0 => (x+h*kx3)),dims=2)) + log.((1 .- (x+h*kx3))./((x+h*kx3) .+ 1))./10
		ka4 = (Δ*exp.(-σ*(dx+h*kd3).^2) .+ δmin - (a+h*ka3))/γ

		dx = (kx1 + 2*kx2 + 2*kx3 + kx4)/6
		da = (ka1 + 2*ka2 + 2*ka3 + ka4)/6

		err = max(maximum(abs.(dx)),maximum(abs.(da)))

		xs = [xs x+h*dx]
		as = [as a+h*da]
	end

	Xs = Matrix{Float64}(undef,n,0)
	As = Matrix{Float64}(undef,n,0)

	for i in 1:c
		Xs = [Xs readdlm("temp/xs_$i.csv",',')]
		rm("temp/xs_$i.csv")
		#As = [As readdlm("temp/as_$i.csv",',')]
		#rm("temp/as_$i.csv")
	end
	Xs = [Xs xs]
	As = [As as]

	return Xs, As, iter, err
end







