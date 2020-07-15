using LinearAlgebra, Statistics, Dates

include("henon-heiles.jl")
include("lorenz.jl")

needed_paths = ["./temp_data/",]
for path in needed_paths
	if !isdir(path)
		mkdir(path)
	end
end

function lyap_hh(xy0::Array{Float64,1}, delta::Float64, n_iter::Int64, n_sample::Int64, h::Float64=.01, lambda::Float64=1.)
	v = 2*rand(2) .- 1
	xy00 = copy(xy0)
	xy00[[1,3]] += delta*v/norm(v)

	c = 0
	cc = 0
	d = Array{Float64,1}()
	XY1 = Array{Float64,2}(undef,4,0)
	XY2 = Array{Float64,2}(undef,4,0)
	while c < n_sample
		c += 1
		
		xys1,E = hh(xy0,false,n_iter,h,lambda)
		xys11,E = hh(xy00,false,n_iter,h,lambda)

		xy1 = xys1[:,end]
		xy11 = xys11[:,end]
		
		if c%1000 == 0
			cc += 1
			writedlm("temp_data/XY1_$cc.csv",XY1,',')
			XY1 = xys1
			writedlm("temp_data/XY2_$cc.csv",XY2,',')
			XY2 = xys11
			writedlm("temp_data/d_$cc.csv",d,',')
			d = [norm(xy1[[1,3]] - xy11[[1,3]]),]
		else
			XY1 = [XY1 xys1]
			XY2 = [XY2 xys11]
			push!(d,norm(xy1[[1,3]] - xy11[[1,3]]))
		end

		xy0 = copy(xy1)
		v = xy11[[1,3]] - xy1[[1,3]]
		v = v/norm(v)
		xy00 = [xy0[1] + delta*v[1],xy11[2],xy0[3]  + delta*v[2],xy11[4]]
	end

	XY11 = Array{Float64,2}(undef,4,0)
	XY22 = Array{Float64,2}(undef,4,0)
	dd = Array{Float64,1}()
	for i in 1:cc
		XY11 = [XY11 readdlm("temp_data/XY1_$i.csv",',')]
		rm("temp_data/XY1_$i.csv")
		XY22 = [XY22 readdlm("temp_data/XY2_$i.csv",',')]
		rm("temp_data/XY2_$i.csv")
		dd = [dd;vec(readdlm("temp_data/d_$i.csv",','))]
		rm("temp_data/d_$i.csv")
	end
	XY11 = [XY11 XY1]
	XY22 = [XY22 XY2]
	dd = [dd;d]

	lyaph = mean(log.(dd./delta)./(h*n_iter))

	return lyaph,dd,XY11,XY22
end


function lyap_lorenz(x0::Array{Float64,1}, delta::Float64, n_iter::Int64, n_sample::Int64, sig::Float64, beta::Float64, rho::Float64, h::Float64)
	v = 2*rand(3) .- 1
	x00 = x0 + delta*v/norm(v)

	c = 0
	cc = 0
	d = Array{Float64,1}()
	X1 = Array{Float64,2}(undef,3,0)
	X2 = Array{Float64,2}(undef,3,0)
	while c < n_sample
		c += 1

		xs = lorenz(x0,n_iter,h,sig,beta,rho)
		x1 = xs[:,end]
		xss = lorenz(x00,n_iter,h,sig,beta,rho)
		x11 = xss[:,end]
		if c%1000 == 0
			cc += 1
			writedlm("temp_data/X1_$cc.csv",X1,',')
			X1 = xs
			writedlm("temp_data/X2_$cc.csv",X2,',')
			X2 = xss
			writedlm("temp_data/d_$cc.csv",d,',')
			d = [norm(x1 - x11),]
		else
			X1 = [X1 xs]
			X2 = [X2 xss]
			push!(d,norm(x1 - x11))
		end

		x0 = copy(x1)
		v = x11 - x1
		v = v/norm(v)
		x00 = x0 + delta*v
	end

	XX1 = Array{Float64,2}(undef,3,0)
	XX2 = Array{Float64,2}(undef,3,0)
	dd = Array{Float64,1}()
	for i in 1:cc
		XX1 = [XX1 readdlm("temp_data/X1_$i.csv",',')]
		rm("temp_data/X1_$i.csv")
		XX2 = [XX2 readdlm("temp_data/X2_$i.csv",',')]
		rm("temp_data/X2_$i.csv")
		dd = [dd;vec(readdlm("temp_data/d_$i.csv",','))]
		rm("temp_data/d_$i.csv")
	end
	XX1 = [XX1 X1]
	XX2 = [XX2 X2]
	dd = [dd;d]

	lyaph = mean(log.(dd./delta)./(h*n_iter))

	return lyaph,dd,XX1,XX2
end










