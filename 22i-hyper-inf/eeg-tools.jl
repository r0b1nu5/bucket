using EDF, DelimitedFiles

function read_eeg(file::String)
	xxx = EDF.read(file)

	s2signal = Dict{String,Vector{Float64}}()

	for sig in xxx.signals
		if typeof(sig) == EDF.Signal{Int16}
			s2signal[sig.header.label] = EDF.decode(sig)
		end
	end

	return s2signal
end

function average_over_zones(s2signal::Dict{String,Vector{Float64}}, s2z::Dict{String,Int64})
	l = maximum(values(s2z))
	T = length(s2signal[collect(keys(s2signal))[1]])

	c = zeros(l)
	asig = zeros(l,T)
	for k in keys(s2signal)
		z = s2z[k]
		c[z] += 1
		asig[z,:] = asig[z,:]*(c[z]-1)/c[z] + s2signal[k]/c[z]
	end

	return asig
end





