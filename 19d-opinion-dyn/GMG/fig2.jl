using DelimitedFiles, PyPlot

X = readdlm("add_data/fig2-X.csv",',')
Y = readdlm("add_data/fig2-Y.csv",',')
Z = readdlm("add_data/fig2-Z.csv",',')

figure("fig2")

PyPlot.contourf(X,Y,Z,LinRange(.4,1.4,50),cmap="YlOrRd")
xlabel("epsilon")
ylabel("Delta")
colorbar(label="effort [%]",ticks=.4:.1:1.4)



X1 = readdlm("add_data/fig2asup-X.csv",',')
Y1 = readdlm("add_data/fig2asup-Y.csv",',')
Z1 = readdlm("add_data/fig2asup-Z.csv",',')
X2 = readdlm("add_data/fig2bsup-X.csv",',')
Y2 = readdlm("add_data/fig2bsup-Y.csv",',')
Z2 = readdlm("add_data/fig2bsup-Z.csv",',')
X3 = readdlm("add_data/fig2csup-X.csv",',')
Y3 = readdlm("add_data/fig2csup-Y.csv",',')
Z3 = readdlm("add_data/fig2csup-Z.csv",',')

figure("fig2sup",(20,4.5))

subplot(1,3,1)
PyPlot.contourf(X1,Y1,Z1,LinRange(0,70,50),cmap="YlOrRd")
xlabel("epsilon")
ylabel("Delta")
colorbar(label="effort [%]",ticks=(0:10:70))

subplot(1,3,2)
PyPlot.contourf(X2,Y2,Z2,LinRange(0,70,50),cmap="YlOrRd")
xlabel("epsilon")
ylabel("Delta")
colorbar(label="effort [%]",ticks=(0:10:70))

subplot(1,3,3)
PyPlot.contourf(X3,Y3,Z3,LinRange(0,.5,50),cmap="YlOrRd")
xlabel("epsilon")
ylabel("Delta")
colorbar(label="effort [%]",ticks=(0:.05:.5))


