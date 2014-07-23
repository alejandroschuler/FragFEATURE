threeLetterAA = ["ALA","CYS","ASP","GLU","PHE","GLY","HIS","ILE","LYS","LEU",
                 "MET","ASN","PRO","GLN","ARG","SER","THR","VAL","TRP","TYR"]

# Define Amino Acid 1-letter to 3-letter conversion
aa1to3dict = {"A":"ALA","C":"CYS","D":"ASP","E":"GLU","F":"PHE","G":"GLY",
              "H":"HIS","I":"ILE","K":"LYS","L":"LEU","M":"MET","N":"ASN",
              "P":"PRO","Q":"GLN","R":"ARG","S":"SER","T":"THR","V":"VAL",
              "W":"TRP","Y":"TYR"}
aa1to3dict["-"] = "GAP"		# Special definition for Gap (-)
aa1to3dict["X"] = "UNK"		# Special definition for Unknown Amino Acid
aa1to3dict["B"] = "ASX"		# Special definition for ASP/ASN
aa1to3dict["Z"] = "GLX"		# Special definition for GLU/GLN

# Define Amin Acid 3-letter to 1-letter conversion
aa3to1dict = {"ALA":"A","CYS":"C","ASP":"D","GLU":"E","PHE":"F","GLY":"G",
              "HIS":"H","ILE":"I","LYS":"K","LEU":"L","MET":"M","ASN":"N",
              "PRO":"P","GLN":"Q","ARG":"R","SER":"S","THR":"T","VAL":"V",
              "TRP":"W","TYR":"Y"}
aa3to1dict["ASX"] = "B"
aa3to1dict["GLX"] = "Z"

aa3to1dict["ACE"] = "X"		# Special definition for ACE
aa3to1dict["NH2"] = "X"		# Special definition for NH2
aa3to1dict["MDO"] = "X"
aa3to1dict["CRO"] = "X"
aa3to1dict["CSW"] = "C"		# Special definition for CYS
aa3to1dict["CSO"] = "C"		# Special definition for CYS
aa3to1dict["MSE"] = "M"		# Special definition for MET
aa3to1dict["CXM"] = "M"		# Special definition for MET
aa3to1dict["FME"] = "M"		# Special definition for MET
aa3to1dict["MLY"] = "K"		# Special definition for LYS
aa3to1dict["MLZ"] = "K"		# Special definition for LYS
aa3to1dict["ELY"] = "K"         # Special definition for LYS
aa3to1dict["NLE"] = "L"		# Special definition for LEU
aa3to1dict["PCA"] = "E"		# Special definition for GLU
aa3to1dict["BME"] = "X"		# Special definition for BME, this is a pdb file error
