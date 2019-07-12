using PyPlot, DelimitedFiles, FFTW

T = 50000
dt = .004
f = 1.
n = 120
fignum = 1

Xs = readdlm("data/uk45_forced_$(f)_$(T)_$(dt).csv",',')

for i in 1:120
	fX = real(fft(Xs[n+i,:]))
	figure(fignum)
	PyPlot.plot((0:25000)./(dt*T),abs.(fX[1:25001]))
end

 #=
for i in [45,]
	fX = real(fft(Xs[n+i,:]))
	
	figure(fignum)
	
	PyPlot.plot((0:25000)./(dt*T),abs.(fX[1:25001])./T,"-b",label="$i")
end


for i in [46,49,66]
	fX = real(fft(Xs[n+i,:]))
	
	figure(fignum)
	
	PyPlot.plot((0:25000)./(dt*T),abs.(fX[1:25001])./T,"-r",label="$i")
end

for i in [47,48,50,54,65,70,117,118]
	fX = real(fft(Xs[n+i,:]))
	
	figure(fignum)
	
	PyPlot.plot((0:25000)./(dt*T),abs.(fX[1:25001])./T,"-g",label="$i")
end

for i in [43,]
	fX = real(fft(Xs[n+i,:]))
	
	figure(fignum)
	
	PyPlot.plot((0:25000)./(dt*T),abs.(fX[1:25001])./T,"-c",label="$i")
end

for i in [42,44]
	fX = real(fft(Xs[n+i,:]))
	
	figure(fignum)
	
	PyPlot.plot((0:25000)./(dt*T),abs.(fX[1:25001])./T,"-y",label="$i")
end

legend()
axis([0,7,1e-4,200])

semilogy()

# =#
