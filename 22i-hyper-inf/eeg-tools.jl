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

function load_and_save_eeg_avg(subjects::Vector{String}, states::Vector{String}, nz::Int64)
	s = readdlm("eeg-data/sensors-$nz.csv",',',String)
	z = readdlm("eeg-data/zones-$nz.csv",',',Int64)
	s2z = Dict{String,Int64}(s[i] => z[i] for i in 1:length(s))
	for su in subjects
		for st in states
			@info "Working on S"*su*"R"*st
			s2signal = read_eeg("eeg-data/S"*su*"R"*st*".edf")
			asig = average_over_zones(s2signal,s2z)
			writedlm("eeg-data/S"*su*"R"*st*"-X.csv",asig,',')
		end
	end

	return nothing
end

function list_all_subjects(n::Int64)
	subjects = String[]
	for i in 1:min(9,n)
		push!(subjects,"00"*string(i))
	end
	for i in 10:min(99,n)
		push!(subjects,"0"*string(i))
	end
	for i in 100:min(999,n)
		push!(subjects,string(i))
	end

	return subjects
end



