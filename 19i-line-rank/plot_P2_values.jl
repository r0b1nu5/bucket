using DelimitedFiles

Asp = readdlm("uk_data/uk_adj_mat.csv",',') .+ 1

n = Int(maximum(Asp))
m = Int(size(Asp)[1]/2)

line_list = Array{Tuple{Int64,Int64},1}()
for i in 1:m
	push!(line_list,(Int64(Asp[2*i-1,1]),Int64(Asp[2*i-1,2])))
end

coord = readdlm("uk_data/uk_grid_coord.csv",',')
bord = readdlm("uk_data/uk_border_coord.csv",',')

figure()
PyPlot.plot(vec(bord[:,1]),vec(bord[:,2]),"-k",linewidth=1)

for i in 1:165
	l = [line_list[i][1],line_list[i][2]]
	p = P2s[i]/losses[i]
#	if 3.1 < p < 3.15
	if p < 6.26
		PyPlot.plot(coord[l,1],coord[l,2],"-y")
#	elseif 3.15 < p < 3.2
	elseif 6.26 < p < 6.34
		PyPlot.plot(coord[l,1],coord[l,2],"-c")
#	elseif 3.2 < p < 3.23
	elseif 6.34 < p < 6.45
		PyPlot.plot(coord[l,1],coord[l,2],"-m")
#	elseif 3.23 < p < 3.28
	elseif 6.45 < p < 6.56
		PyPlot.plot(coord[l,1],coord[l,2],"-b")
#	elseif 3.28 < p < 3.33
	elseif 6.56 < p < 6.64
		PyPlot.plot(coord[l,1],coord[l,2],"-r")
#	elseif 3.33 < p < 3.38
	elseif 6.64 < p < 6.75
		PyPlot.plot(coord[l,1],coord[l,2],"-g")
#	elseif 3.38 < p < 4
	elseif 6.75 < p < 7
		PyPlot.plot(coord[l,1],coord[l,2],"-c")
	else
		PyPlot.plot(coord[l,1],coord[l,2],":k")
	end
end

PyPlot.plot(vec(coord[:,1]),vec(coord[:,2]),"ok")

