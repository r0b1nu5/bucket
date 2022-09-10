using DelimitedFiles

js = ["nature",
	"pnas",
	"science",
	"lancet",#"lancet_letters",
	"neng_j_med",#"neng_j_med_letters",
	"plant_cell",
	"j_phys_chem_a",
	"ieee_trans_autom_control",
	"energy",
	"chaos",#"chaos_91-99","chaos_00-08","chaos_09-17",
	"siam",
	"ann_math",
	"prl",
	"prd"]

ds = ["pl","plc","yule"]

ns = 5000

for j in js
	@info "=============================================="
	@info j
	for d in ds
		dat = readdlm("analysis/"*j*"_"*d*"_params_$(ns).csv",',')
		p = readdlm("analysis/"*j*"_"*d*"_p-gof_$(ns).csv",',')

		@info d*": params = ($(round.(dat[1:end],digits=3))), p-val = $(round(p[1],digits=4))"
	end
end



