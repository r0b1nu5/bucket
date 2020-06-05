# Computes the number of representative for each party.
#
# np: Number of votes for party +1.
# nn: Number of votes for party -1.
# m: Number of representatives for the state.

function result(x::Array{Float64,1}, m::Int64)
	n = length(x)

	np = sum(x .>0)
	nn = n - np

	pp = np/(np+nn)
	pn = nn/(np+nn)

	mp = round(Int64,m*pp)
	mn = round(Int64,m*pn)

	while mp + mn > m
		@info "There's a tie, one representative is attributed at random."
		if rand() > .5
			mn -= 1
		else
			mp -= 1
		end
	end
	while mp + mn < m
		@info "There's a tie, one representative is attributed at random."
		if rand() > .5
			mn += 1
		else
			mp += 1
		end
		
	end

	return mp, mn
end


# Computes the margin before each party looses a representative

function result_margin(x::Array{Float64,1}, m::Int64)
	n = length(x)

	np = sum(x .> 0)
	nn = n - np

	pp = np/n
	pn = nn/n

	mp = round(Int64,m*pp)
	mn = round(Int64,m*pn)

	while mp + mn > m
		@info "There's a tie, one representative is attributed at random."
		if rand() > .5
			mn -= 1
		else
			mp -= 1
		end
	end
	while mp + mn < m
		@info "There's a tie, one representative is attributed at random."
		if rand() > .5
			mn += 1
		else
			mp += 1
		end
		
	end

	imp = floor(m*pp)
	if imp == mp
		imp -= 1
	end
	imn = floor(m*pn)
	if imn == mn
		imn -= 1
	end

	inp = (2*imp+1)*n/(2*m)
	inn = (2*imn+1)*n/(2*m)

	margp = np - inp
	margn = nn - inn

	return mp, margp, mn, margn
end






