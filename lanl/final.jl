using DelimitedFiles, PyPlot, FFTW, Statistics

# Estimates the forcings frequency, based on the Fourier Transform of the times series of the phase frequencies. A peak is identified by dividing the value of each Fourier mode by the median value of its 2*Df neighbors.
#
# INPUT:
# Xs: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
# dt: Time step.
# Df: Radius of data over which the median in computed. 
#
# OUTPUT:
# fh: Estimate of the forcing's frequency.

function get_fh_fourier(Xs::Array{Float64,2}, dt::Float64, Df::Int=10)
	nn,T = size(Xs)
	n = Int(nn/2)

	fX = zeros(Complex{Float64},nn,T)
	for i in 1:n
		fX[i,:] = fft(Xs[n+i,:]).*dt./pi
	end

	nfX = norm.(fX)
	
	mfX = zeros(n,T-2*Df)
	for i in 1:n
		for j in Df+1:T-Df
			mfX[i,j-Df] = nfX[i,j]/median([nfX[i,j-Df:j-1];nfX[i,j+1:j+Df]])
		end
	end

	fs = (0:T-1)./(dt*T)
	df = fs[2]-fs[1]

	freqs = Array{Float64,1}()
	maxs = Array{Float64,1}()

	for i in 1:n
		ma,id = findmax(mfX[i,1:Int(T/2)])

		ma2,id2 = findmax([mfX[i,1:id-1];0.;mfX[i,id+1:Int(T/2)]])

		if (ma2 > .5*ma) && (abs(id - id2) == 1)
			push!(freqs,mean([fs[Df+id],fs[Df+id2]]))
		elseif (ma2 > .8*ma)
			push!(freqs,NaN)
			@info "Fourier Transform: inconclusive!"
		else
			push!(freqs,fs[Df+id])
		end
	end

	fh = median(freqs)
	tf = sum((abs.(freqs .- fh) .> 2*df))

	if tf/n > .05
		@info "Fourier Transform: no clear result."

		return NaN
	else
		return fh
	end
end


# Estimates the forcing's frequency based on the autocorrelation of the times series of the phase angles.
#
# INPUT:
# Xs: Time series of the phase angles (rows 1:n) and of the phase frequencies (n+1:2*n).
# dt: Time step.
# p: Proportion of the measure time over which the autocorrelation is computed.
#
# OUTPUT: 
# fh: Estimate of the forcing's frequency.

function get_fh_autocorr(Xs::Array{Float64,2}, dt::Float64, p::Float64=.1)
	nn,T = size(Xs)
	n = Int(nn/2)

	ac = Array{Float64,2}(undef,n,0)
	leng = floor(Int,p*T)
	for i in 1:leng
		XX = Xs[1:n,1:end-i].*Xs[1:n,i+1:end]
		ac = [ac vec(sum(XX,dims=2))./(T-i)]
	end
	@info "Autocorrelation computed."

	fh = Array{Float64,1}()
	for i in 1:n
		tops = Array{Int64,1}()
		bots = Array{Int64,1}()

		for j in 2:leng-1
			if (ac[i,j-1] > ac[i,j]) && (ac[i,j+1] > ac[i,j])
				push!(bots,j)
			elseif (ac[i,j-1] < ac[i,j]) && (ac[i,j+1] < ac[i,j])
				push!(tops,j)
			end
		end
		dtops = tops[2:end] - tops[1:end-1]
		dbots = bots[2:end] - bots[1:end-1]
		del = median([dtops;dbots])

		push!(fs,1/(del*dt))
	end

	return median(fh)
end










