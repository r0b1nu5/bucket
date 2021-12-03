using PyPlot, DelimitedFiles, Statistics

name_cmap = "viridis"
cmap = get_cmap(name_cmap)
# #=
col1 = cmap(1/6)
col2 = cmap(3/6)
col3 = cmap(5/6)

run = 6179
n = 18

t = (0:n)*2π/n

x = Array(LinRange(-1,1,50))
z1 = Array(LinRange(0,40,50))
z2 = Array(LinRange(0,4π,50))
X = ones(50)*x'
Y = sqrt.(1 .- X.^2)
Z1 = z1*ones(1,50)
Z2 = z2*ones(1,50)

θ0 = vec(readdlm("temp_data/th0f_$(run).csv",','))[[4:18;1:3]]
d0 = mod.(θ0[[2:n;1]] - θ0 .+ π,2π) .- π
θ0 = [(θ0[1] + sum(d0[1:i])) for i in 0:n]
θ2 = vec(readdlm("temp_data/th2f_$(run).csv",','))[[4:18;1:3]]
d2 = mod.(θ2[[2:n;1]] - θ2 .+ π,2π) .- π
θ2 = [(θ2[1] + sum(d2[1:i])) for i in 0:n]

figure()
subplot(2,1,1)
PyPlot.plot([0,20],[-π,-π],"--k")
PyPlot.plot([0,20],[0,0],"--k")
PyPlot.plot([0,20],[π,π],"--k")
#PyPlot.plot(1:n+1,θ0 .- mean(θ0) .- 2π,"-o",color=col2)
PyPlot.plot(n+1:2*n+1,θ0 .- mean(θ0),"-o",color=col2)
PyPlot.plot(1:n+1,θ0 .- mean(θ0),"-o",color=col2)
PyPlot.plot(1-n:1,θ0 .- mean(θ0),"-o",color=col2)
#PyPlot.plot(1:n+1,θ0 .- mean(θ0) .+ 2π,"-o",color=col2)
axis([0,n+2,-π-.5,π+.5])
xticks(1:2:n+1,[1:2:n;1])
yticks(π*[-1,-.5,0,.5,1])

subplot(2,1,2)
PyPlot.plot([0,20],[-π,-π],"--k")
PyPlot.plot([0,20],[0,0],"--k")
PyPlot.plot([0,20],[π,π],"--k")
PyPlot.plot(1-n:1,θ2 .- mean(θ2) .+ 1.1*π/2 .- 2π,"-o",color=col1)
PyPlot.plot(1:n+1,θ2 .- mean(θ2) .+ 1.1*π/2 .- 2π,"-o",color=col1)
PyPlot.plot(n+1:2*n+1,θ2 .- mean(θ2) .+ 1.1*π/2 .- 2π,"-o",color=col1)
PyPlot.plot(1-n:1,θ2 .- mean(θ2) .+ 1.1*π/2,"-o",color=col1)
PyPlot.plot(1:n+1,θ2 .- mean(θ2) .+ 1.1*π/2,"-o",color=col1)
PyPlot.plot(n+1:2*n+1,θ2 .- mean(θ2) .+ 1.1*π/2,"-o",color=col1)
PyPlot.plot(1-n:1,θ2 .- mean(θ2) .+ 1.1*π/2 .+ 2π,"-o",color=col1)
PyPlot.plot(1:n+1,θ2 .- mean(θ2) .+ 1.1*π/2 .+ 2π,"-o",color=col1)
PyPlot.plot(n+1:2*n+1,θ2 .- mean(θ2) .+ 1.1*π/2 .+ 2π,"-o",color=col1)
axis([0,n+2,-π-.5,π+.5])
xticks(1:2:n+1,[1:2:n;1])
yticks(π*[-1,-.5,0,.5,1])

#=
figure()

PyPlot.plot_surface(.95*X,.95*Y,Z1,color="grey",alpha=.2,shade=false)
PyPlot.plot_surface(.95*X,-.95*Y,Z1,color="grey",alpha=.2,shade=false)

PyPlot.plot3D(cos.(θ0),sin.(θ0),Array(1:n+1),"-o",color=col2)
PyPlot.plot3D(cos.(θ2),sin.(θ2),Array(1:n+1) .+ 20,"-o",color=col1)

figure()

PyPlot.plot_surface(.95*X,.95*Y,Z2,color="grey",alpha=.2,shade=false)
PyPlot.plot_surface(.95*X,-.95*Y,Z2,color="grey",alpha=.2,shade=false)

PyPlot.plot3D(cos.(t),sin.(t),θ0,"-o",color=col2)
PyPlot.plot3D(cos.(t),sin.(t),θ2 .+ 2π,"-o",color=col1)

=#


run = 6786
n = 18

t = (0:n)*2π/n

θ0 = vec(readdlm("temp_data/th0f_$(run).csv",','))[[3:18;1:2]]
d0 = mod.(θ0[[2:n;1]] - θ0 .+ π,2π) .- π
θ0 = [(θ0[1] + sum(d0[1:i])) for i in 0:n]
θ2 = vec(readdlm("temp_data/th1f_$(run).csv",','))[[3:18;1:2]]
d2 = mod.(θ2[[2:n;1]] - θ2 .+ π,2π) .- π
θ2 = [(θ2[1] + sum(d2[1:i])) for i in 0:n]

figure()
subplot(2,1,1)
PyPlot.plot([0,20],[-π,-π],"--k")
PyPlot.plot([0,20],[0,0],"--k")
PyPlot.plot([0,20],[π,π],"--k")
#PyPlot.plot(1:n+1,θ0 .- mean(θ0) .- 2π,"-o",color=col2)
PyPlot.plot(n+1:2*n+1,θ0 .- mean(θ0),"-o",color=col2)
PyPlot.plot(1:n+1,θ0 .- mean(θ0),"-o",color=col2)
PyPlot.plot(1-n:1,θ0 .- mean(θ0),"-o",color=col2)
#PyPlot.plot(1:n+1,θ0 .- mean(θ0) .+ 2π,"-o",color=col2)
axis([0,n+2,-π-.5,π+.5])
xticks(1:2:n+1,[1:2:n;1])
yticks(π*[-1,-.5,0,.5,1])

subplot(2,1,2)
PyPlot.plot([0,20],[-π,-π],"--k")
PyPlot.plot([0,20],[0,0],"--k")
PyPlot.plot([0,20],[π,π],"--k")
PyPlot.plot(1-n:1,θ2 .- mean(θ2) .- 2π,"-o",color=col1)
PyPlot.plot(1:n+1,θ2 .- mean(θ2) .- 2π,"-o",color=col1)
PyPlot.plot(n+1:2*n+1,θ2 .- mean(θ2) .- 2π,"-o",color=col1)
PyPlot.plot(1-n:1,θ2 .- mean(θ2),"-o",color=col1)
PyPlot.plot(1:n+1,θ2 .- mean(θ2),"-o",color=col1)
PyPlot.plot(n+1:2*n+1,θ2 .- mean(θ2),"-o",color=col1)
PyPlot.plot(1-n:1,θ2 .- mean(θ2) .+ 2π,"-o",color=col1)
PyPlot.plot(1:n+1,θ2 .- mean(θ2) .+ 2π,"-o",color=col1)
PyPlot.plot(n+1:2*n+1,θ2 .- mean(θ2) .+ 2π,"-o",color=col1)
axis([0,n+2,-π-.5,π+.5])
xticks(1:2:n+1,[1:2:n;1])
yticks(π*[-1,-.5,0,.5,1])


