using Graphs

function my_er(n::Int64, p::Float64, tol::Float64=1e-6, max_try::Int64=100)
	l2 = 0.
	c = 0

	while l2 < tol && c < max_try
		c += 1
		g = erdos_renyi(n,p)
		global B = incidence_matrix(g; oriented=true)
		l2 = laplacian_spectrum(g)[2]
	end

	return Float64.(B)
end


# TBC
#


