using PyPlot, LinearAlgebra

function jacobian(th::Array{Float64,1}, LL::Array{Float64,2}=zeros(1,1))
	n = length(th)

	if size(LL) == (1,1)
		L = n*diagm(0 => ones(n)) - ones(n,n)
		A = 1 .- diagm(0 => ones(n))
	else
		L = LL
		A = -L.*(1 .- diagm(0 => ones(n)))
	end

	J0 = A.*cos.(th*ones(1,n) - ones(n)*th')
	J = J0 - diagm(0 => vec(sum(J0,dims=2)))

	return J
end

function get_quartiles(x::Matrix{Float64})
	n,m = size(x)

	q0 = [quantile(x[:,i],0.) for i in 1:m]
	q1 = [quantile(x[:,i],.25) for i in 1:m]
	q2 = [quantile(x[:,i],.5) for i in 1:m]
	q3 = [quantile(x[:,i],.75) for i in 1:m]
	q4 = [quantile(x[:,i],1.) for i in 1:m]

	return q0,q1,q2,q3,q4
end

function plot_quart(qs::Tuple{Any,Any,Any,Any,Any}, τs::Vector{Float64}, col::String="C0")
	q0,q1,q2,q3,q4 = qs

	n = length(τs)

	PyPlot.plot(1 ./τs,q0,"x",color=col)
	PyPlot.plot(1 ./τs,q4,"x",color=col)
	for i in 1:n
		PyPlot.plot([1/τs[i],1/τs[i]],[q1[i],q3[i]],"-",color=col,linewidth=2.)
	end
	PyPlot.plot(1 ./τs,q2,"o",color=col)

	return nothing
end



