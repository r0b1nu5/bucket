using Distributions

function bi_gauss(n::Int64, m1::Float64, s1::Float64, m2::Float64, s2::Float64)
	n1 = floor(Int,n/2)
	n2 = n - n1

	y1 = rand(Normal(m1,s1),n1)
	y2 = rand(Normal(m2,s2),n2)

	return [y1;y2]
end

function gauss(n::Int64, m::Float64, s::Float64)
	return rand(Normal(m,s),n)
end




