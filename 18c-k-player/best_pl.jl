include("ml_zipf.jl")
include("journals.jl")
include("my_histo.jl")
include("gof.jl")

using DelimitedFiles, PyPlot

js = ["prl"]

for j in js
	num = Array{Float64,1}(vec(readdlm("./data/"*j*".txt",'\t')[2:end-2]))
	mi = minimum(num)
	ma = maximum(num)
	mas = intersect(num)
	mat = copy(mas)
	ss = Array{Float64,1}()
	go = Array{Float64,1}()
	
	while mat[end] > mi+9
		tes = (1:length(num)) .* (num .<= mat[end])
		m = maximum(tes)
		push!(ss,ml_clauset(num[1:m]))
		push!(go,gof(num[1:m]))

		mat = mat[1:end-1]
	end
	
	figure()
	PyPlot.plot((mi+10):mas[end],ss)
	figure()
	PyPlot.plot((mi+10):mas[end],go)
end

(max_go,ind) = findmax(go)
max_s = ss[ind]



