using PyPlot

α = .5

x1 = LinRange(-1.5*π,1.5*π,200)
x2 = LinRange(-π/2-α,π/2-α,100)
x3 = LinRange(-π/2+α,π/2+α,100)

y1 = sin.(x1)
y2 = sin.(x1 .- α) .+ sin(α)
y3 = sin.(-x1 .- α) .+ sin(α)
y4 = sin.(x3 .- α) .+ sin(α)
y5 = sin.(-x2 .- α) .+ sin(α)

PyPlot.plot([-1.5*π,1.5*π],[0,0],"k")
PyPlot.plot([0,0],[-1.2,1.5],"k")

PyPlot.plot(x1,y2,color="C1")
PyPlot.plot(x1,y3,color="C2")
PyPlot.plot(x3,y4,color="C1",linewidth=4.)
PyPlot.plot(x2,y5,color="C2",linewidth=4.)
PyPlot.plot(x1,y1,"--",color="C0")


axis([-1.5*π,1.5*π,-1.2,1.6])
