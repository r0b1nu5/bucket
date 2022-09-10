using PyPlot

include("../get_rgb.jl")

n_run = 5
n_srun = 1000
ntw = "ntw3"
n = 3
ee = -.00005

its = Array{Int64,1}()
th0s = Array{Float64,2}(undef,n,0)

for i in 1:n_run
	global its,th0s
	itt = Int.(readdlm("data/"*ntw*"_$(i)_its.csv",','))
	th0t = readdlm("data/"*ntw*"_$(i)_th0s.csv",',')
	its = [its;itt]
	th0s = [th0s th0t]
end

mi,imi = findmin(its)
ma,ima = findmax(its)
rgb = get_rgb(mi,ma)

figure("pd",(6,6))

ids = setdiff(Array(1:length(its)),[imi,ima])

for i in ids
	PyPlot.plot(th0s[1,i]-th0s[3,i],th0s[2,i]-th0s[3,i],".",color=rgb[its[i]])
end
PyPlot.plot(th0s[1,imi]-th0s[3,imi],th0s[2,imi]-th0s[3,imi],".",color=rgb[its[imi]],label="$(mi)")
PyPlot.plot(th0s[1,ima]-th0s[3,ima],th0s[2,ima]-th0s[3,ima],".",color=rgb[its[ima]],label="$(ma)")



title("Phase diagram, ε = $(ee)")
xlabel("θ1 - θ3")
ylabel("θ2 - θ3")
legend()

