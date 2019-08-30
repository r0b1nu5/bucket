using PyPlot,Statistics

include("kuramoto_q.jl")
include("../default_color_cycle.jl")

n_iter = 10000
q = 3

ini = [3*pi*rand(Float64,2,n_iter).-3*pi/2;zeros(1,n_iter)]

fin1 = Array{Float64,2}(undef,3,0)
finq = Array{Float64,2}(undef,3,0)

for i in 1:n_iter
	global ini,fin1
	
	th = vec(kuramoto_q(1,zeros(3),1.,ini[:,i],false))
	
	fin1 = [fin1 th.-th[end]]
end

mix1 = minimum(round.(Int,fin1[1,:]))
max1 = maximum(round.(Int,fin1[1,:]))
miy1 = minimum(round.(Int,fin1[2,:]))
may1 = maximum(round.(Int,fin1[2,:]))

c = 0
figure(1)

for x in mix1:max1
	for y in miy1:may1
		t1 = round.(Int,fin1[1,:]./(2*pi)) .== x
		t2 = round.(Int,fin1[2,:]./(2*pi)) .== y
		idx = setdiff(t1.*t2.*(1:n_iter),[0,])
		figure(1)
		PyPlot.plot(ini[1,idx],ini[2,idx],"x")
	end
	global fin1,n_iter,ini
end

for i in 1:n_iter
	global ini,finq

	th = vec(kuramoto_q(q,zeros(3),1/q,ini[:,i],false))
	
	finq = [finq th.-th[end]]
end

mix2 = minimum(round.(Int,finq[1,:]))
max2 = maximum(round.(Int,finq[1,:]))
miy2 = minimum(round.(Int,finq[2,:]))
may2 = maximum(round.(Int,finq[2,:]))

c = 0
figure(2)

for x in mix2:max2
	for y in miy2:may2
		t1 = round.(Int,finq[1,:]*q/(2*pi)) .== x
		t2 = round.(Int,finq[2,:]*q/(2*pi)) .== y
		idx = setdiff(t1.*t2.*(1:n_iter),[0,])
		if length(idx) > 0
			c += 1
		end
		figure(2)
		PyPlot.plot(ini[1,idx],ini[2,idx],"x",color=def_col[1+c%9])
	end
	global finq,n_iter,ini,c
end




