using PyPlot

include("tools.jl")
include("ksakaguchi.jl")

n = 3

B = [1 1 0;-1 0 1;0 -1 -1.]
C = [1 -1 1.]
L = B*B'

ω = ones(n)
ϕ = .3
T = 500

figure(1,(14,5))

figure(2,(14,5))

col1 = "C0"
col2 = "C1"
col3 = "C2"

x0 = 2π*rand(n) .- π

xx = ksakaguchi(L,ω,x0,ϕ,true,false,.01,-1.,T)
yy = ksakaguchi(L,ω,x0 + 2π*[1,0,0],ϕ,true,false,.01,-1.,T)
zz = ksakaguchi(L,ω,x0 .+ 1.,ϕ,true,false,.01,-1.,T)

figure(1)
PyPlot.plot3D(xx[1][1,:],xx[1][2,:],xx[1][3,:],"-",color=col1)
PyPlot.plot3D(yy[1][1,:],yy[1][2,:],yy[1][3,:],"--",color=col1)
PyPlot.plot3D(zz[1][1,:],zz[1][2,:],zz[1][3,:],":",color=col1)

figure(2)
subplot(3,1,1)
PyPlot.plot((0:T)*.01,xx[1][1,:],color=col1)
PyPlot.plot((0:T)*.01,yy[1][1,:],color=col1)
PyPlot.plot((0:T)*.01,zz[1][1,:],color=col1)
subplot(3,1,2)
PyPlot.plot((0:T)*.01,xx[1][2,:],color=col1)
PyPlot.plot((0:T)*.01,yy[1][2,:],color=col1)
PyPlot.plot((0:T)*.01,zz[1][2,:],color=col1)
subplot(3,1,3)
PyPlot.plot((0:T)*.01,xx[1][3,:],color=col1)
PyPlot.plot((0:T)*.01,yy[1][3,:],color=col1)
PyPlot.plot((0:T)*.01,zz[1][3,:],color=col1)

x1 = 2π*rand(n) .- π
ω = rand(n)

xx = ksakaguchi(L,ω,x1,ϕ,true,false,.01,-1.,T)
yy = ksakaguchi(L,ω,x1 + 2π*[1,0,0],ϕ,true,false,.01,-1.,T)
zz = ksakaguchi(L,ω,x1 .+ 1.,ϕ,true,false,.01,-1.,T)

figure(1)
PyPlot.plot3D(xx[1][1,:],xx[1][2,:],xx[1][3,:],"-",color=col2)
PyPlot.plot3D(yy[1][1,:],yy[1][2,:],yy[1][3,:],"--",color=col2)
PyPlot.plot3D(zz[1][1,:],zz[1][2,:],zz[1][3,:],":",color=col2)

figure(2)
subplot(3,1,1)
PyPlot.plot((0:T)*.01,xx[1][1,:],color=col2)
PyPlot.plot((0:T)*.01,yy[1][1,:],color=col2)
PyPlot.plot((0:T)*.01,zz[1][1,:],color=col2)
subplot(3,1,2)
PyPlot.plot((0:T)*.01,xx[1][2,:],color=col2)
PyPlot.plot((0:T)*.01,yy[1][2,:],color=col2)
PyPlot.plot((0:T)*.01,zz[1][2,:],color=col2)
subplot(3,1,3)
PyPlot.plot((0:T)*.01,xx[1][3,:],color=col2)
PyPlot.plot((0:T)*.01,yy[1][3,:],color=col2)
PyPlot.plot((0:T)*.01,zz[1][3,:],color=col2)

x2 = 2π*rand(n) .- π
ω = 10*rand(3)

xx = ksakaguchi(L,ω,x2,ϕ,true,false,.01,-1.,T)
yy = ksakaguchi(L,ω,x2 + 2π*[1,0,0],ϕ,true,false,.01,-1.,T)
zz = ksakaguchi(L,ω,x2 .+ 1.,ϕ,true,false,.01,-1.,T)

figure(1)
PyPlot.plot3D(xx[1][1,:],xx[1][2,:],xx[1][3,:],"-",color=col3)
PyPlot.plot3D(yy[1][1,:],yy[1][2,:],yy[1][3,:],"--",color=col3)
PyPlot.plot3D(zz[1][1,:],zz[1][2,:],zz[1][3,:],":",color=col3)

figure(2)
subplot(3,1,1)
PyPlot.plot((0:T)*.01,xx[1][1,:],color=col3)
PyPlot.plot((0:T)*.01,yy[1][1,:],color=col3)
PyPlot.plot((0:T)*.01,zz[1][1,:],color=col3)
subplot(3,1,2)
PyPlot.plot((0:T)*.01,xx[1][2,:],color=col3)
PyPlot.plot((0:T)*.01,yy[1][2,:],color=col3)
PyPlot.plot((0:T)*.01,zz[1][2,:],color=col3)
subplot(3,1,3)
PyPlot.plot((0:T)*.01,xx[1][3,:],color=col3)
PyPlot.plot((0:T)*.01,yy[1][3,:],color=col3)
PyPlot.plot((0:T)*.01,zz[1][3,:],color=col3)

figure(1)
xlabel("x1")
ylabel("x2")
zlabel("x3")
