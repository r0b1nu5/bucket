using PyPlot

include("tools.jl")

n = 3

ω = zeros(n)
ω = rand(n)
#ω = [1.,1.,1.]
ϕ = .3

B = [1 1 0;-1 0 1;0 -1 -1.]
C = [1 -1 1.]
L = B*B'

cols = [(161,207,239)./255,(76,163,224)./255,(31,119,180)./255]

#=
xres = 1000
xmin = -4π
kmin = floor(Int64,(xmin+π)/(2π))
xmax = 6π
kmax = ceil(Int64,(xmax-π)/(2π))
yres = 500
ymin = -2π
lmin = floor(Int64,(ymin+π)/(2π))
ymax = 5π/2
lmax = ceil(Int64,(ymax-π)/(2π))


xs = repeat(Vector(LinRange(xmin,xmax,xres))',yres,1)
ys = repeat(Vector(LinRange(ymin,ymax,yres)),1,xres)

u = zeros(yres,xres)
v = zeros(yres,xres)
#q = zeros(Int64,yres,xres)

for i in 1:yres
	for j in 1:xres
		f = (B*(sin.(B'*[xs[j],ys[i],0.] .+ ϕ) .- sin(ϕ)))[1:2]
		u[i,j] = f[1]
		v[i,j] = f[2]
#		q[i,j] = winding([xs[j],ys[i],0.],[1,2,3])
	end
end

speed = sqrt.(u.^2 + v.^2) 
uu = u.*(speed .> .05)
vv = v.*(speed .> .05)
=#


res = 200
xs = repeat((Vector(LinRange(-π,π,res+1))[1:res])',res,1)
ys = repeat(Vector(LinRange(-π,π,res+1))[1:res],1,res)
ks = -1:2
ls = -1:2
Xs = zeros(length(ls)*res,length(ks)*res)
Ys = zeros(length(ls)*res,length(ks)*res)
xmin = -2π+1
xmax = 4π 
ymin = -2π+2.5
ymax = 7π/2
for k in 1:length(ks)
	for l in 1:length(ls)
		Xs[(l-1)*res+1:l*res,(k-1)*res+1:k*res] = xs[1:res,1:res] .+ 2π*ks[k]
		Ys[(l-1)*res+1:l*res,(k-1)*res+1:k*res] = ys[1:res,1:res] .+ 2π*ls[l]
	end
end

l1,l2 = size(Xs)

u = zeros(l1,l2)
v = zeros(l1,l2)

for i in 1:l1
	for j in 1:l2
		f = (ω - B*(sin.(B'*[Xs[i,j],Ys[i,j],0.] .- ϕ) .+ sin(ϕ)))[1:2]
		u[i,j] = f[1]
		v[i,j] = f[2]
	end
end

speed = sqrt.(u.^2 + v.^2) 


x0 = [-π,-π,0.,π,π,0.,-π]
y0 = [-π,0.,π,π,0.,-π,-π]
xm1 = [0.,π,π,0.]
ym1 = [-π,0.,-π,-π]
xp1 = [-π,-π,0.,-π]
yp1 = [0.,π,π,0.]

figure("flow",(14,5))

#=
PyPlot.plot([xmin,xmax],[π,π],"--",color="C7",linewidth=.5)
PyPlot.plot([xmin,xmax],[-π,-π],"--",color="C7",linewidth=.5)
PyPlot.plot([-3π,-3π],[ymin,ymax],"--",color="C7",linewidth=.5)
PyPlot.plot([-π,-π],[ymin,ymax],"--",color="C7",linewidth=.5)
PyPlot.plot([π,π],[ymin,ymax],"--",color="C7",linewidth=.5)
PyPlot.plot([3π,3π],[ymin,ymax],"--",color="C7",linewidth=.5)
PyPlot.plot([5π,5π],[ymin,ymax],"--",color="C7",linewidth=.5)
=#

for k in ks #kmin:kmax
	for l in ls #lmin:lmax
		PyPlot.fill(x0.+2π*k,y0.+2π*l,color=cols[1])
		PyPlot.fill(xm1.+2π*k,ym1.+2π*l,color=cols[2])
		PyPlot.fill(xp1.+2π*k,yp1.+2π*l,color=cols[3])
	end
end

ρ = .9

#PyPlot.streamplot(xs,ys,u,v,density=4,linewidth=.5,color="black",arrowsize=.7)#,maxlength=.5)
PyPlot.streamplot(Xs,Ys,u,v,density=4.5,linewidth=2*(ρ*speed./maximum(speed).+(1-ρ)),color="black",arrowsize=.8)#,maxlength=.5)
#PyPlot.streamplot(xs,ys,u,v,density=4,linewidth=.8,color=speed,cmap="plasma",arrowsize=.8,maxlength=.2,integration_direction="forward")

axis([xmin,xmax,ymin,ymax])
xlabel("x1")
ylabel("x2")
