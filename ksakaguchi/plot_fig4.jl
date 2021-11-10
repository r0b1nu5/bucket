using PyPlot, DelimitedFiles

name_cmap = "viridis"
cmap = get_cmap(name_cmap)
col1 = cmap(1/6)
col2 = cmap(3/6)
col3 = cmap(5/6)

# #=
run = 6179

n = 18

α = .4

ω = [1.;zeros(5);-1.;zeros(11)]

 #=
include("gen_fig4.jl")
# =#

# #=
ps = vec(readdlm("temp_data/fig4_ps_$(run).csv",','))

λ0 = readdlm("temp_data/fig4_l0_$(run).csv",',')
λ1 = readdlm("temp_data/fig4_l1_$(run).csv",',')
λ2 = readdlm("temp_data/fig4_l2_$(run).csv",',')

ω0 = vec(readdlm("temp_data/fig4_w0_$(run).csv",','))
ω1 = vec(readdlm("temp_data/fig4_w1_$(run).csv",','))
ω2 = vec(readdlm("temp_data/fig4_w2_$(run).csv",','))

# =#

figure()

# #=
subplot(1,2,1)
PyPlot.plot(ps,λ1[end-1,:],color=col3,label="sol. 3")
PyPlot.plot(ps,λ0[end-1,:],color=col2,label="sol. 2")
PyPlot.plot(ps,λ2[end-1,:],color=col1,label="sol. 1")
xlabel("p")
ylabel("λ_2")
legend()

subplot(1,2,2)
# =#

 #=
for i in 1:n
	subplot(1,4,1)
	PyPlot.plot(ps,λ0[i,:])
	subplot(1,4,2)
	PyPlot.plot(ps,λ1[i,:])
	subplot(1,4,3)
	PyPlot.plot(ps,λ2[i,:])
end

subplot(1,4,4)
# =#

PyPlot.plot(ps,-ω1,color=col3)
PyPlot.plot(ps,-ω0,color=col2)
PyPlot.plot(ps,-ω2,color=col1)
xlabel("p")
ylabel("ω_s")

# =#

# #=
run = 6786

n = 18

α = .4

ω = [1.;zeros(8);-1.;zeros(8)]

 #=
include("gen_fig4.jl")
# =#

# #=
ps = vec(readdlm("temp_data/fig4_ps_$(run).csv",','))

λ0 = readdlm("temp_data/fig4_l0_$(run).csv",',')
λ1 = readdlm("temp_data/fig4_l1_$(run).csv",',')
λ2 = readdlm("temp_data/fig4_l2_$(run).csv",',')

ω0 = vec(readdlm("temp_data/fig4_w0_$(run).csv",','))
ω1 = vec(readdlm("temp_data/fig4_w1_$(run).csv",','))
ω2 = vec(readdlm("temp_data/fig4_w2_$(run).csv",','))

# =#

figure()

# #=
subplot(1,2,1)
#PyPlot.plot(ps,λ2[end-1,:],color=col3,label="sol. 3")
PyPlot.plot(ps,λ0[end-1,:],color=col2,label="sol. 2")
PyPlot.plot(ps,λ1[end-1,:],color=col1,label="sol. 1")
xlabel("p")
ylabel("λ_2")
legend()

subplot(1,2,2)
# =#

 #=
for i in 1:n
	subplot(1,4,1)
	PyPlot.plot(ps,λ0[i,:])
	subplot(1,4,2)
	PyPlot.plot(ps,λ1[i,:])
	subplot(1,4,3)
	PyPlot.plot(ps,λ2[i,:])
end

subplot(1,4,4)
# =#

#PyPlot.plot(ps,-ω2,color=col3)
PyPlot.plot(ps,-ω0,color=col2)
PyPlot.plot(ps,-ω1,color=col1)
xlabel("p")
ylabel("ω_s")

# =#
