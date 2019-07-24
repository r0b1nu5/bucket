using LinearAlgebra

# s: number of nodes on the side of the triangle

function build_tri(s::Int64, plot::Bool = false)
	v1 = [cos(-pi/6),sin(-pi/6)]
	v2 = [0., 1.]
	
	n = Int.(s*(s+1)/2)
	
	A = zeros(n,n)
	
	xy = Array{Float64,2}(undef,0,2)
	
	c = 0
	
	for i in 1:s
		for j in 1:i
			c += 1
			
			xy = [xy; ([0., 0.] + (i-1)*v1 + (j-1)*v2)']
			
			if j > 1
				A[c,c-1] = A[c-1,c] = 1.
				
				if i > 1
					A[c,c-i] = A[c-i,c] = 1.
				end
			end
			
			if i > 1 && j < i
				A[c,c-i+1] = A[c-i+1,c] = 1.
			end
		end
	end
	
	L = diagm(0 => vec(sum(A,dims=2))) - A
	
	if plot
		figure("Triangle $n")
		for i in 1:n-1
			for j in i+1:n
				if L[i,j] != 0
					PyPlot.plot(xy[[i,j],1],xy[[i,j],2],"-ok")
				end
			end
			PyPlot.text(xy[i,1],xy[i,2],"$i")
		end
		PyPlot.text(xy[n,1],xy[n,2],"$n")
	end
	
	return L, xy
end
			
			




