using Distributions

function big_rand(n1::Int64, m1::Float64, s1::Float64, n2::Int64, m2::Float64, s2::Float64)

	x0 = sort([rand(Normal(m1,s1),n1); rand(Normal(m2,s2),n2)])

	return x0
end


