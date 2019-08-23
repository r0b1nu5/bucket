# List of journals considered.
# For each journal, we removed the Anonymous author, as well as the authors whose number of articles had not all data due to limits in file size on WoS.

# #=
journals_short = ["nature",
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
#
	"prl",#"prl_58-77","prl_78-97","prl_98-17",
# =#
#journals_short = ["prd",]
#
# #=
	"bmj",
	"food_chem_tox",
	"medical_ped_onc",
	"pre",
	"molecular_bio_evo",
	"scientometrics",
	"j_mat_chem_a",
	"bull_ams"]
# =#
used_journals = ["chaos","ieee_trans_autom_control","lancet","nature","neng_j_med","pnas","prd","pre","prl","science","bmj","medical_ped_onc","scientometrics"]

journals_full = Dict{String,String}()

journals_full["ann_math"] = "Annals of Mathematics"
journals_full["bmj"] = "BMJ"
journals_full["bull_ams"] = "Bulletin of the American Mathematical Society"
journals_full["chaos"] = "Chaos"
journals_full["chaos_91-99"] = "Chaos (1991-99)"
journals_full["chaos_95-03"] = "Chaos (1995-2003)"
journals_full["chaos_00-08"] = "Chaos (2000-08)"
journals_full["chaos_04-12"] = "Chaos (2004-12)"
journals_full["chaos_09-17"] = "Chaos (2009-17)"
journals_full["energy"] = "Energy"
journals_full["food_chem_tox"] = "Food and Chemical Toxicology"
journals_full["ieee_trans_autom_control"] = "IEEE Transactions on Automatic Control"
journals_full["j_phys_chem_a"] = "Journal of Physical Chemistry A"
journals_full["j_mat_chem_a"] = "Journal of Material Chemistry A"
journals_full["lancet"] = "Lancet"
journals_full["lancet_letters"] = "Letters of Lancet"
journals_full["medical_ped_onc"] = "Medical and Pediatric Oncology"
journals_full["molecular_bio_evo"] = "Molecular Biology and Evolution"
journals_full["nature"] = "Nature"
journals_full["nature_comm"] = "Nature Communications"
journals_full["neng_j_med"] = "New England Journal of Medicine"
journals_full["neng_j_med_letters"] = "Letters of New England Journal of Medicine"
journals_full["plant_cell"] = "Plant Cell"
journals_full["pnas"] = "PNAS"
journals_full["prd"] = "Physical Review D"
journals_full["pre"] = "Physical Review E"
journals_full["prl"] = "Physical Review Letters"
journals_full["prl_58-77"] = "Physical Review Letters (1958-77)"
journals_full["prl_68-87"] = "Physical Review Letters (1968-87)"
journals_full["prl_78-97"] = "Physical Review Letters (1978-97)"
journals_full["prl_88-07"] = "Physical Review Letters (1988-2007)"
journals_full["prl_98-17"] = "Physical Review Letters (1998-2017)"
journals_full["science"] = "Science"
journals_full["scientometrics"] = "Scientometrics"
journals_full["siam"] = "SIAM Journal on Applied Mathematics"

journals_code = Dict{String,String}()

journals_code["chaos"] = "A"
journals_code["energy"] = "E-bis"
journals_code["ieee_trans_autom_control"] = "B"
journals_code["lancet"] = "C"
journals_code["nature"] = "D"
journals_code["neng_j_med"] = "E"
journals_code["pnas"] = "F"
journals_code["prd"] = "G-bis"
journals_code["pre"] = "G"
journals_code["prl"] = "H"
journals_code["science"] = "I"
journals_code["bmj"] = "J"
journals_code["medical_ped_onc"] = "K"
journals_code["scientometrics"] = "L"

journals_code["food_chem_tox"] = "M"

journals_colors = Dict{String,Tuple{Array{Float64,1},Array{Float64,1}}}()

journals_colors["nature"] = ([1.,0.,0.],[0.,0.,0.])
journals_colors["pnas"] = ([1.,.5,0.],[0.,0.,0.])
journals_colors["science"] = ([1.,1.,0.],[0.,0.,0.])
journals_colors["lancet"] = ([.5,1.,0.],[0.,0.,0.])
journals_colors["neng_j_med"] = ([0.,1.,0.],[0.,0.,0.])
journals_colors["plant_cell"] = ([0.,1.,.5],[0.,0.,0.])
journals_colors["j_phys_chem_a"] = ([0.,1.,1.],[0.,0.,0.])
journals_colors["ieee_trans_autom_control"] = ([0.,.5,1.],[0.,0.,0.])
journals_colors["energy"] = ([0.,0.,1.],[0.,0.,0.])
journals_colors["chaos"] = ([.5,0.,1.],[0.,0.,0.])
journals_colors["siam"] = ([1.,0.,1.],[0.,0.,0.])
journals_colors["ann_math"] = ([1.,0.,.5],[0.,0.,0.])

journals_colors["prl"] = ([.3,.3,.3],[0.,0.,0.])
journals_colors["prd"] = ([.6,.6,.6],[0.,0.,0.])
# + cross-figure

## Appendix...
journals_colors["bmj"] = ([1.,0.,0.],[0.,0.,0.])
journals_colors["food_chem_tox"] = ([1.,.5,0.],[0.,0.,0.])
journals_colors["medical_ped_onc"] = ([1.,1.,0.],[0.,0.,0.])
journals_colors["pre"] = ([.5,1.,0.],[0.,0.,0.])
journals_colors["molecular_bio_evo"] = ([0.,1.,0.],[0.,0.,0.])
journals_colors["scientometrics"] = ([0.,1.,.5],[0.,0.,0.])
journals_colors["j_mat_chem_a"] = ([0.,1.,1.],[0.,0.,0.])
journals_colors["bull_ams"] = ([0.,.5,1.],[0.,0.,0.])




