# List of journals considered.
# For each journal, we removed the Anonymous author, as well as the authors whose number of articles had not all data due to limits in file size on WoS.

# #=
journals_short = ["nature",
	"pnas",
	"science",
	"lancet",#"lancet_letters",
	"neng_j_med",#"neng_j_med_letters",
	"plant_cell",
	"j_acs",
	"ieee_trans_autom_control",
	"energy",
	"chaos",#"chaos_91-99","chaos_00-08","chaos_09-17",
	"siam",
	"ann_math",
	"prl",
	"prd"]
#
#	"prl",#"prl_58-77","prl_78-97","prl_98-17",
# =#
#journals_short = ["prl",]
#journals_short = ["prd",]
#
 #=
	"prd",
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

journals_red_short = [
		      "nature_1950",
		      "pnas_1950",
		      "science_1940",
		      "lancet_1910",
		      "neng_j_med_1950",
		      "plant_cell_2000",
		      "j_acs_1930",
		      "ieee_trans_autom_control_2000",
		      "energy_2005",
		      "chaos",
		      "siam",
		      "ann_math"
		      ]


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

journals_code["ann_math"] = "AMA"
journals_code["bmj"] = "BMJ"
journals_code["bull_ams"] = "BAM"
journals_code["chaos"] = "CHA"
journals_code["energy"] = "ENE"
journals_code["energy_2005"] = "ENE*"
journals_code["food_chem_tox"] = "FCT"
journals_code["ieee_trans_autom_control"] = "TAC"
journals_code["ieee_trans_autom_control_2000"] = "TAC*"
journals_code["j_phys_chem_a"] = "PCA"
journals_code["j_mat_chem_a"] = "MCA"
journals_code["j_acs"] = "ACS"
journals_code["j_acs_1930"] = "ACS*"
journals_code["lancet"] = "LAN"
journals_code["lancet_1910"] = "LAN*"
journals_code["medical_ped_onc"] = "MPO"
journals_code["molecular_bio_evo"] = "MBE"
journals_code["nature"] = "NAT"
journals_code["nature_1950"] = "NAT*"
journals_code["nature_comm"] = "NAC"
journals_code["neng_j_med"] = "NEM"
journals_code["neng_j_med_1950"] = "NEM*"
journals_code["plant_cell"] = "PLC"
journals_code["plant_cell_2000"] = "PLC*"
journals_code["pnas"] = "PNA"
journals_code["pnas_1950"] = "PNA*"
journals_code["prd"] = "PRD"
journals_code["pre"] = "PRE"
journals_code["prl"] = "PRL"
journals_code["science"] = "SCI"
journals_code["science_1940"] = "SCI*"
journals_code["scientometrics"] = "SCM"
journals_code["siam"] = "SIA"

### Colors from the "category20" panel.
journals_colors = Dict{String,Tuple{Array{Float64,1},Array{Float64,1}}}()

journals_colors["nature"] = (parse.(Int64,["1f","77","b4"],base = 16)./255,[0.,0.,0.])
journals_colors["nature_1950"] = (parse.(Int64,["1f","77","b4"],base = 16)./255,[0.,0.,0.])
journals_colors["pnas"] = (parse.(Int64,["ff","7f","0e"],base = 16)./255,[0.,0.,0.])
journals_colors["pnas_1950"] = (parse.(Int64,["ff","7f","0e"],base = 16)./255,[0.,0.,0.])
journals_colors["science"] = (parse.(Int64,["2c","a0","2c"],base = 16)./255,[0.,0.,0.])
journals_colors["science_1940"] = (parse.(Int64,["2c","a0","2c"],base = 16)./255,[0.,0.,0.])
journals_colors["lancet"] = (parse.(Int64,["1f","77","b4"],base = 16)./255,[0.,0.,0.])
journals_colors["lancet_1910"] = (parse.(Int64,["1f","77","b4"],base = 16)./255,[0.,0.,0.])
journals_colors["neng_j_med"] = (parse.(Int64,["94","67","bd"],base = 16)./255,[0.,0.,0.])
journals_colors["neng_j_med_1950"] = (parse.(Int64,["94","67","bd"],base = 16)./255,[0.,0.,0.])
journals_colors["plant_cell"] = (parse.(Int64,["8c","56","4b"],base = 16)./255,[0.,0.,0.])
journals_colors["plant_cell_2000"] = (parse.(Int64,["8c","56","4b"],base = 16)./255,[0.,0.,0.])
journals_colors["j_phys_chem_a"] = (parse.(Int64,["e3","77","c2"],base = 16)./255,[0.,0.,0.])
journals_colors["j_acs"] = (parse.(Int64,["e3","77","c2"],base = 16)./255,[0.,0.,0.])
journals_colors["j_acs_1930"] = (parse.(Int64,["e3","77","c2"],base = 16)./255,[0.,0.,0.])
journals_colors["ieee_trans_autom_control"] = (parse.(Int64,["7f","7f","7f"],base = 16)./255,[0.,0.,0.])
journals_colors["ieee_trans_autom_control_2000"] = (parse.(Int64,["7f","7f","7f"],base = 16)./255,[0.,0.,0.])
journals_colors["energy"] = (parse.(Int64,["bc","bd","22"],base = 16)./255,[0.,0.,0.])
journals_colors["energy_2005"] = (parse.(Int64,["bc","bd","22"],base = 16)./255,[0.,0.,0.])
journals_colors["chaos"] = (parse.(Int64,["17","be","cf"],base = 16)./255,[0.,0.,0.])
journals_colors["siam"] = (parse.(Int64,["ae","c7","e8"],base = 16)./255,[0.,0.,0.])
journals_colors["ann_math"] = (parse.(Int64,["ff","bb","78"],base = 16)./255,[0.,0.,0.])

journals_colors["prl"] = (parse.(Int64,["98","df","8a"],base = 16)./255,[0.,0.,0.])
journals_colors["prd"] = (parse.(Int64,["ff","98","96"],base = 16)./255,[0.,0.,0.])
# + cross-figure

## Appendix...
journals_colors["bmj"] = (parse.(Int64,["1f","77","b4"],base = 16)./255,[0.,0.,0.])
journals_colors["food_chem_tox"] = (parse.(Int64,["ff","7f","0e"],base = 16)./255,[0.,0.,0.])
journals_colors["medical_ped_onc"] = (parse.(Int64,["2c","a0","2c"],base = 16)./255,[0.,0.,0.])
journals_colors["pre"] = (parse.(Int64,["d6","27","28"],base = 16)./255,[0.,0.,0.])
journals_colors["molecular_bio_evo"] = (parse.(Int64,["94","67","bd"],base = 16)./255,[0.,0.,0.])
journals_colors["scientometrics"] = (parse.(Int64,["8c","56","4b"],base = 16)./255,[0.,0.,0.])
journals_colors["j_mat_chem_a"] = (parse.(Int64,["e3","77","c2"],base = 16)./255,[0.,0.,0.])
journals_colors["bull_ams"] = (parse.(Int64,["7f","7f","7f"],base = 16)./255,[0.,0.,0.])




