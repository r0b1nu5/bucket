using SpecialFunctions

# Computes the volume of the qf-basin, conditioned on the initial winding number to be qi and assuming the distribtion is normal (Wiley et al. 2006).
function vol_qf_normal(qf::Int64, qi::Int64, n::Int64, mfi::Float64, sfi::Float64, mif::Float64, sif::Float64)
	Vi = vol_wcell(qi,n) # To be implemented

	return (erf((qf-mfi*qi+.5)/(sfi*sqrt(2n))) - erf((qf-mfi*qi-.5)/(sfi*sqrt(2*n))))/(erf((qi-mif*qf+.5)/(sif*sqrt(2*n))) - erf((qi-mif*qf-.5)/(sif*sqrt(2*n))))
end

# Computes the volume of the qf-basin, conditioned on the initial winding number being qi and assuming the distribution is exponential (Delabays et al. 2017).
function vol_qf_exp(qf::Int64, qi::Int64, n::Int64, mfi::Float64, lfi::Float64, mif::Float64, lif::Float64)
	Vi = vol_wcell(qi,n)
	qmax = floor(Int64,n/4)

	return (exp(-lfi*abs(qf-mfi*qi)/sqrt(n))*(1-exp(lfi/sqrt(n)))*(2*cosh(lif*mif*qf/sqrt(n))*exp(-lif*qmax/sqrt(n)) - exp(-lif*(mif*qf-ceil(mif*qf))/sqrt(n)) - exp(-lif*(floor(mif*qf)-mif*qf)/sqrt(n))))/(exp(-lif*abs(qi-mif*qf)/sqrt(n))*(1-exp(lif/sqrt(n)))*(2*cosh(lfi*mfi*qi/sqrt(n))*exp(-lfi*qmax/sqrt(n)) - exp(-lfi*(mfi*qi-ceil(mfi*qi))/sqrt(n)) - exp(-lfi*(floor(mfi*qi)-mfi*qi)/sqrt(n))))
end

function vol_qf_exp_old(qf::Int64, qi::Int64, n::Int64, mfi::Float64, lfi::Float64, mif::Float64, lif::Float64)
	Vi = vol_wcell(qi,n)

	return (sinh(lfi/(2*sqrt(n)))*exp(-lfi*abs(qf-mfi*qi)/sqrt(n)))/(sinh(lif/(2*sqrt(n)))*exp(-lif*abs(qi-mif*qf)/sqrt(n)))
end

# Volume of the qi-winding cell in a n-cycle.
function vol_wcell(qi::Int64, n::Int64)
	return 1.
end



