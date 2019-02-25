# List of journals considered.
# For each journal, we removed the Anonymous author, as well as the authors whose number of articles had not all data due to limits in file size on WoS.

journals_short = ["bmj",
	"chaos","chaos_91-99","chaos_00-08","chaos_09-17",
	"food_chem_tox",
	"ieee_trans_autom_control",
	"lancet","lancet_letters",
	"medical_ped_onc",
	"nature",
	"neng_j_med","neng_j_med_letters",
	"pnas",
	"prl","prl_58-77","prl_78-97","prl_98-17",
	"pre",
	"science",
	"scientometrics"]

used_journals = ["chaos","ieee_trans_autom_control","lancet","nature","neng_j_med","pnas","pre","prl","science","bmj","medical_ped_onc","scientometrics"]

journals_full = {
	"bmj": "BMJ",
	"chaos": "Chaos",
	"chaos_91-99": "Chaos (1991-99)",
	"chaos_95-03": "Chaos (1995-2003)",
	"chaos_00-08": "Chaos (2000-08)",
	"chaos_04-12": "Chaos (2004-12)",
	"chaos_09-17": "Chaos (2009-17)",
	"food_chem_tox": "Food and Chemical Toxicology",
	"ieee_trans_autom_control": "IEEE Transactions on Automatic Control",
	"lancet": "Lancet",
	"lancet_letters": "Letters of Lancet",
	"medical_ped_onc": "Medical and Pediatric Oncology",
	"nature": "Nature",
	"neng_j_med": "New England Journal of Medicine",
	"neng_j_med_letters": "Letters of New England Journal of Medicine",
	"pnas": "PNAS",
	"pre": "Physical Review E",
	"prl": "Physical Review Letters",
	"prl_58-77": "Physical Review Letters (1958-77)",
	"prl_68-87": "Physical Review Letters (1968-87)",
	"prl_78-97": "Physical Review Letters (1978-97)",
	"prl_88-07": "Physical Review Letters (1988-2007)",
	"prl_98-17": "Physical Review Letters (1998-2017)",
	"science": "Science",
	"scientometrics": "Scientometrics"
}

journals_code = {
	"chaos": "A",
	"ieee_trans_autom_control": "B",
	"lancet": "C",
	"nature": "D",
	"neng_j_med": "E",
	"pnas": "F",
	"pre": "G",
	"prl": "H",
	"science": "I",
	"bmj": "J",
	"medical_ped_onc": "K",
	"scientometrics": "L",
	"food_chem_tox": "M"
}

journals_colors = {
	"chaos": ([1.,0.,0.],[1.,0.,1.]),
	"ieee_trans_autom_control": ([0.,1.,0.],[0.,.8,0.]),
	"lancet": ([0.,0.,1.],[0.,1.,1.]),
	"nature": ([.9,.9,0.],[1.,.5,0.]),
	"neng_j_med": ([1.,0.,1.],[1.,0.,0.]),
	"pnas": ([0.,1.,1.],[0.,0.,1.]),
	"pre": ([0.,0.,0.],[.5,.5,.5]),
	"prl": ([1.,.5,0.],[.9,.9,0.]),
	"science": ([.5,.5,.5],[0.,0.,0.]),
	"bmj": ([148/255,0.,211/255],[0.,0.,0.]),
	"medical_ped_onc": ([0.,.5,0.],[0.,0.,0.]),
	"scientometrics": ([.25,.87,.81],[0.,0.,0.])
}


