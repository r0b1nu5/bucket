using LinearAlgebra

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


