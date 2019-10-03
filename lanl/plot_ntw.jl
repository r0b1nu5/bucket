using PyPlot, DelimitedFiles

function plot_ntw(ntw::String)
	if ntw == "ntw5"
		n = 5

		x = cos.(2*pi/5*(0:4))
		y = sin.(2*pi/5*(0:4))

		L = readdlm("data/ntw5_lap_mat.csv",',')

		for i in 1:4
			for j in i+1:5
				if abs(L[i,j]) > 0
					PyPlot.plot([x[i],x[j]],[y[i],y[j]],"-k")
				end
			end
		end

		PyPlot.plot(x,y,"ok")

	elseif ntw == "ntw10"
		n = 10

		x = cos.(2*pi/10*(0:9))
		y = sin.(2*pi/10*(0:9))

		L = readdlm("data/ntw10_lap_mat.csv",',')

		for i in 1:9
			for j in i+1:10
				if abs(L[i,j]) > 0
					PyPlot.plot([x[i],x[j]],[y[i],y[j]],"-k")
				end
			end
		end

		PyPlot.plot(x,y,"ok")

	elseif ntw == "uk"
		n = 120

		xy = readdlm("uk/uk_grid_coord.csv",',')
		x = xy[:,1]
		y = xy[:,2]

		bord = readdlm("uk/uk_border_coord.csv",',')
		xb = bord[:,1]
		yb = bord[:,2]
		
		L = readdlm("uk/uk_lap_mat.csv",',')

		PyPlot.plot([xb;xb[1]],[yb;yb[1]],"-k",linewidth=.5)

		for i in 1:n-1
			for j in i+1:n
				if abs(L[i,j]) > 0
					PyPlot.plot([x[i],x[j]],[y[i],y[j]],"-k",linewidth=1.)
				end
			end
		end

		PyPlot.plot(x,y,"ok",markersize=3.)

		axis([-.1,1.1,0,1])
	end

	xticks([])
	yticks([])
end



