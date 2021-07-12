using PyPlot

function get_ex(x,d,hgm,hgp)
	f2s = LinRange(hgm,hgp-x,200)

	Y = Array{Float64,1}()
	for f2 in f2s
		f1 = f2 + x
		y1 = hi(f1)
		y2 = hi(f2)

		push!(Y,y1-y2)
	end

	figure()
	PyPlot.plot(f2s,Y)
	PyPlot.plot([f2s[1],f2s[end]],[d*x,d*x])
end


function test_solv(x,D,hgm,hgp)
	n = length(x)

	t = Array{Bool,1}()
	
	for i in 1:n
		s = sign(x[i])

		if s == 1.
			y2 = hi(hgm)
			y1 = hi(hgm + x[i])
			dy1 = y1 - y2

			y3 = hi(x[i])
			dy2 = y3

			if (dy1 - D[i,i]*x[i])*(dy2 - D[i,i]*x[i]) < 0.
				push!(t,true)
			else
				push!(t,false)
			end
		else
			y2 = hi(hgm - x[i])
			y1 = hi(hgm)
			dy1 = y1 - y2

			y3 = hi(-x[i])
			dy2 = -y3

			if (dy1 - D[i,i]*x[i])*(dy2 - D[i,i]*x[i]) < 0.
				push!(t,true)
			else
				push!(t,false)
			end
		end
	end

	return t
end






