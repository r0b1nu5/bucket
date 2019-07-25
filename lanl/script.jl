ti = time()

fX45 = fft(Xs[120+45,:])
fX114 = fft(Xs[120+114,:])

figure()
subplot(211)
PyPlot.plot((0:T-1)./(dt*T),real.(fX45))
subplot(212)
PyPlot.plot((0:T-1)./(dt*T),real.(fX114))

Ah, fh = find_A_n_f(Xs,dt)

ah, ph = locate_f_n_lag(Xs,dt,fh,Ah)

@info "Total time: $(time() - ti)"




