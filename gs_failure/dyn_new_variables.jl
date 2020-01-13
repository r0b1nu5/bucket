using PyPlot


(x0,y0) = 2pi*rand(2)

err = 1000.
h = .0001
b = .91

x1 = copy(x0)
y1 = copy(y0)

xs = Array{Float64,1}()
ys = Array{Float64,1}()

iter = 0

while err > 1e-8 && iter < 100000
	global iter,xs,ys,x1,y1,b,h,err
	iter += 1
	
	push!(xs,x1)
	push!(ys,y1)

	kx1 = -3*b - 3*sin(x1)*cos(y1)
	ky1 = -cos(x1)*sin(y1)
	kx2 = -3*b - 3*sin(x1+h/2*kx1)*cos(y1+h/2*ky1)
	ky2 = -cos(x1+h/2*kx1)*sin(y1+h/2*ky1)
	kx3 = -3*b - 3*sin(x1+h/2*kx2)*cos(y1+h/2*ky2)
	ky3 = -cos(x1+h/2*kx2)*sin(y1+h/2*ky2)
	kx4 = -3*b - 3*sin(x1+h*kx3)*cos(y1+h*ky3)
	ky4 = -cos(x1+h*kx3)*sin(y1+h*ky3)

	dx = (kx1+2*kx2+2*kx3+kx4)/6
	dy = (ky1+2*ky2+2*ky3+ky4)/6

	x1 += h*dx
	y1 += h*dy

	err = max(abs(dx),abs(dy))
end

xs = mod.(xs,2pi)
ys = mod.(ys,2pi)

dx = xs[2:end] - xs[1:end-1]
dy = ys[2:end] - ys[1:end-1]

xjumps = setdiff((1:length(xs)-1).*(abs.(dx) .> 1.),[0.,])
yjumps = setdiff((1:length(ys)-1).*(abs.(dx) .> 1.),[0.,])
jumps = [0;union(xjumps,yjumps);length(xs)]

for k in 2:length(jumps)
	PyPlot.plot(xs[(jumps[k-1]+1):(jumps[k])],ys[(jumps[k-1]+1):(jumps[k])],color=col)
end


