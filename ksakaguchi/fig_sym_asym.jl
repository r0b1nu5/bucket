using PyPlot

α = .1
γm = -π/2 + α
γp = π/2 - α

δ = .7

η = .5 # slope of the afine extention h_η

x = LinRange(γm,γp,300)


figure()

subplot(1,3,1)

x1 = x[280]
y1 = cos(α)*sin(x1)
s1 = cos(α)*cos(x1)
dx1 = [-δ,δ]
dy1 = s1*dx1
PyPlot.plot(x1 .+ dx1,y1 .+ dy1,"--k",linewidth=1.)

x2 = x[20]
y2 = cos(α)*sin(x2)
s2 = cos(α)*cos(x2)
dx2 = [-δ,δ]
dy2 = s2*dx2
PyPlot.plot(x2 .+ dx2,y2 .+ dy2,"--k",linewidth=1.)

PyPlot.plot(x,cos(α)*sin.(x),color="C0")
PyPlot.plot([γp,γp + δ],[0,η*δ] .+ cos(α)*sin(γp),"--",color="C0")
PyPlot.plot([γm,γm - δ],[0,-η*δ] .+ cos(α)*sin(γm),"--",color="C0")


subplot(1,3,2)

PyPlot.plot(x,sin.(x .- α) .+ sin(α),color="C1")
PyPlot.plot([γp,γp + δ],[0,η*δ] .+ sin(γp-α) .+ sin(α),"--",color="C1")
PyPlot.plot([γm,γm - δ],[0,-η*δ] .+ sin(γm-α) .+ sin(α),"--",color="C1")

PyPlot.plot(x,sin.(-x .- α) .+ sin(α), color="C2")
PyPlot.plot([γp,γp + δ],[0,-η*δ] .+ sin(-γp-α) .+ sin(α),"--",color="C2")
PyPlot.plot([γm,γm - δ],[0,η*δ] .+ sin(-γm-α) .+ sin(α),"--",color="C2")


subplot(1,3,3)

x3 = x[280]
y3 = sin(α)*(1 - cos(x3))
s3 = sin(α)*sin(x3)
dx3 = [-δ,δ]
dy3 = s3*dx3
PyPlot.plot(x3 .+ dx3,y3 .+ dy3,"--k",linewidth=1.)

x4 = x[20]
y4 = sin(α)*(1 - cos(x4))
s4 = sin(α)*sin(x4)
dx4 = [-δ,δ]
dy4 = s4*dx4
PyPlot.plot(x4 .+ dx4,y4 .+ dy4,"--k",linewidth=1.)

PyPlot.plot(x,sin(α)*(1 .- cos.(x)),color="C3")
PyPlot.plot([γp,γp + δ],[1,1]*sin(α)*(1-cos(γp)),"--",color="C3")
PyPlot.plot([γm,γm - δ],[1,1]*sin(α)*(1-cos(γm)),"--",color="C3")



