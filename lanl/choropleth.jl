#=

NOTE: You must install this data set for this gist: https://www2.census.gov/geo/tiger/TIGER2017/COUNTY/tl_2017_us_county.zip

Copyright 2021 Jacob Zelko (aka TheCedarPrince)
Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
=#

using Plots, RecipesBase, Shapefile

state_shape_file = "data_utk/states_shapes/tl_2021_us_state.shp"
states = Shapefile.Table(state_shape_file)

#=

02 - Alaska
15 - Hawaii
60 - American Samoa
66 - Guam
69 - Northern Mariana Islands
72 - Puerto Rico
78 - Virgin Islands

=#

skip = ["02", "15", "60", "66", "69", "72", "78"] # For creating pretty output

function selectShapes(table, skipped)
    geoms = empty(Shapefile.shapes(table))
    for row in table
        if !(row.STATEFP in skipped)
            push!(geoms, Shapefile.shape(row))
        end
    end
    return geoms
end

ss = selectShapes(states, skip)

states_xy = Vector{Tuple{Vector{Float64},Vector{Float64}}}()
for s in ss
	x = Vector{Float64}()
	y = Vector{Float64}()
	for p in s.points
		push!(x,p.x)
		push!(y,p.y)
	end
	push!(states_xy,(x,y))
end

for s in states_xy
	PyPlot.plot(s[1],s[2],"-k",linewidth=.5)
end



#|> x -> plot(x, axis = ([], false))
