function get_rgb(qmin::Int,qmax::Int)
	rgb = Dict{Int,Tuple{Float64,Float64,Float64}}()
	
	for q in qmin:qmax
		col = abs(q)/qmax
		if col <= .25
			rgb[q] = (0,col*4,1)
		elseif col <= .5
		        rgb[q] = (0,1,2-col*4)
		elseif col <= .75
		        rgb[q] = (col*4-2,1,0)
		else
		        rgb[q] = (1,4-col*4,0)
		end
	end
	
	return rgb
	
end

# gives a color scale between vmin and vmax


function get_rgb(v::Float64,vmin::Float64,vmax::Float64)
	
	col = (v-vmin)/(vmax-vmin)
	
	if col <= .25
		return (0,col*4,1)
	elseif col <=.5
		return (0,1,2-col*4)
	elseif col <=.75
		return (col*4-2,1,0)
	else
		return (1,4-col*4,0)
	end

end
