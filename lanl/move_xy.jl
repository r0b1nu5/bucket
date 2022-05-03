using LinearAlgebra

# Tries to move coordinates so they do not overlap...

function move_xy(x0,y0,λ,δ=1e-4,T=100,tol=1e-4)
	n = length(x0)
	x = x0
	y = y0

	err = 1000.
	c = 0
	
	while err > tol && c < T
		c += 1
		@info "c = $c"

		ndx = zeros(n,n)
		dx = zeros(n,n)
		dy = zeros(n,n)

		for i in 1:n
			for j in 1:n
				ndx[i,j] = norm(x[i]-x[j])
				dx[i,j] = x[i] - x[j]
				dy[i,j] = y[i] - y[j]
			end
		end
		
		xd = (exp.(-λ.*ndx).*sign.(dx))*ones(n)
		yd = (exp.(-λ.*ndx).*sign.(dy))*ones(n)

		x += δ.*xd
		y += δ.*yd

		err = max(maximum(abs.(xd)),maximum(abs.(yd)))
	end

	return x,y
end




