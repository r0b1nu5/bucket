function lv(x0::Vector{Float64}, r::Vector{Float64}, θ::Vector{Float64}, B::Matrix{Float64}, h::Float64=.01, save::Bool=false, tol::Float64=1e-6, maxiter::Int64=10000)
	x = copy(x0)
	xs = copy(x0)

	corr = 1e10
	iter = 0

	while corr > tol && iter < maxiter
		iter += 1
		if iter%100 == 0
			@info "iter: $iter, log-err: $(round(log(corr))))"
		end

		k1 = f_lv(x,r,θ,B)
		k2 = f_lv(x + h*k1/2,r,θ,B)
		k3 = f_lv(x + h*k2/2,r,θ,B)
		k4 = f_lv(x + h*k3,r,θ,B)

		dx = (k1+2*k2+2*k3+k4)/6
		corr = maximum(abs.(dx))

		x += h*dx
		if save 
			xs = [xs x]
		end
	end

	if save 
		return xs
	else
		return x
	end
end


function f_lv(x::Vector{Float64}, r::Vector{Float64}, θ::Vector{Float64}, B::Matrix{Float64})
	return x.*(r - θ.*x + B*x)
end






