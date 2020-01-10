using DelimitedFiles


function load_ntw(ntw::String, ee::Float64)
	if ntw == "ntw3"
		L = [1. 0. -1.;0. 1. -1.;-1. -1. 2.]
		P = [1-ee;1-ee;2*ee-2]
	end

	return L,P
end




