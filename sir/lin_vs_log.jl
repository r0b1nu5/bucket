using Statistics, DelimitedFiles, PyPlot

I_ref = vec(readdlm("data/I_ref_0229.csv",','))[2:end]

x = Array(1:length(I_ref))

di = 14

d0 = Date(2020,2,29)

r_lin = Array{Float64,1}()
r_log = Array{Float64,1}()
dates = Array{String,1}()

for i in 1:length(I_ref)-di
	global d0
	d0 += Dates.Day(1)
	dd = d0 + Dates.Day(di)
	push!(dates,"$(Dates.day(d0)).$(Dates.month(d0)) - $(Dates.day(dd)).$(Dates.month(dd))")
	push!(r_lin,cor(x[i:i+di],I_ref[i:i+di]))
	push!(r_log,cor(x[i:i+di],log.(I_ref[i:i+di])))
end

figure()
subplot(2,1,1)
PyPlot.plot(x[1:length(I_ref)-di],r_lin,"-o",label="Linear")
PyPlot.plot(x[1:length(I_ref)-di],r_log,"-o",label="Exponential")

xticks(1:length(I_ref)-di,dates,rotation=70)
ylabel("Pearson correlation")

legend()


