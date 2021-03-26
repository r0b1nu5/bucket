using PyPlot, DelimitedFiles

run = 0

αs = vec(readdlm("temp_data/αs_$(run).csv",','))
β0min = vec(readdlm("temp_data/β0min_$(run).csv",','))
β0max = vec(readdlm("temp_data/β0max_$(run).csv",','))
β1min = vec(readdlm("temp_data/β1min_$(run).csv",','))
β1max = vec(readdlm("temp_data/β1max_$(run).csv",','))

PyPlot.plot(αs,β0min,color="C0")
PyPlot.plot(αs,β0max,color="C0")
PyPlot.plot(αs,β1min,color="C1")
PyPlot.plot(αs,β1max,color="C1")
title("Run: $run")
xlabel("α")
ylabel("width")


