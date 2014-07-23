import string
import re
###############################################################################
# FF_line object contains all information stored in a line of a FF file
class FF_line():
    def __init__(self,line):
        line = line.strip().split("\t")
        self.ID = line[0]	  # example: Env_1a2l.1_0
        self.props = line[1:-6]	  # numeric properties of microenvironment
        self.x = float(line[-5])  # microenvironment x coordinate
        self.y = float(line[-4])  # microenvironment y coordinate 
        self.z = float(line[-3])  # microenvironment z coordinate
        self.sep = line[-2]	  # separator string before annotation info
        self.info = line[-1] 	  # sring of user specified info
        self.entry = self.ID.split("_")[1]  # example: 1a2l.1

    def update_info(self,new_info):
        self.info = new_info

    def __str__(self):
        string = "%s\t%s\t%s\t%.3f\t%.3f\t%.3f\t%s\t%s" %(self.ID, 
                 "\t".join(self.props), self.sep, self.x, self.y, self.z,
                 self.sep, self.info)
	return string


###############################################################################
# FF_point object contains all non-property information stored in a line of a 
# FF file
class FF_point():
    def __init__(self,line):
        line = line.strip().split("\t")
        self.ID = line[0]
        self.x = float(line[1])
        self.y = float(line[2])
        self.z = float(line[3])
        self.info = line[4]

	self.aaInfo = self.info.split("|")[0]
	self.aaName = self.aaInfo.split(":")[0]
	self.atmName = self.aaInfo.split(":")[1]
	try:
	    self.aaNum = int(self.aaInfo.split(":")[2])
	    self.iCode = " "
	except:
	    self.aaNum = int(self.aaInfo.split(":")[2][:-1])
	    self.iCode = self.aaInfo.split(":")[2][-1]
	self.chainID = self.aaInfo.split(":")[3]

        self.file = self.ID.split("_")[1]
        self.tag  = self.file.split(".")[0].upper() + "." + self.chainID

	self.liginfo = self.info.split("|")[2:-1]
	self.ligIDs = [ligatom.split(":")[0] for ligatom in self.liginfo]
	self.ligIDs = list(set(self.ligIDs))
	self.ligtags = ["%s.%s" %(ligatom.split(":")[0], ligatom.split(":")[2])
                        for ligatom in self.liginfo]

    def __str__(self):
        return "%s\t%.3f\t%.3f\t%.3f\t%s" %(self.ID, self.x, self.y, self.z, 
               self.info)

    def distance(self,self2):
        dx = pow(self.x - self2.x, 2)
        dy = pow(self.y - self2.y, 2)
        dz = pow(self.z - self2.z, 2)
        return pow(dx + dy + dz, 0.5)

    def update_aaName(self,new_aaName):
	self.aaName = new_aaName
	self.aaInfo = ":".join(map(str, [self.aaName, self.atmName, self.aaNum,
                                         self.chainID]))
	self.info = self.info.split("|")
	self.info = self.aaInfo + "|" + "|".join(self.info[1:])

    def update_atmName(self,new_atmName):
	self.atmName = new_atmName
	self.aaInfo = ':'.join(map(str,[self.aaName, self.atmName, self.aaNum, 
                                        self.chainID]))
	self.info = self.info.split("|")
	self.info = self.aaInfo + "|" + "|".join(self.info[1:])

    def update_info(self,new_info):
        self.info = new_info


###############################################################################
# PDB_line object contains all information stored in a line of a pdb file
# Reference: http://www.wwpdb.org/documentation/format32/v3.2.html
class PDB_line:
    def __init__(self, line):                           # Output formatting
        line = line.strip()                                
        self.lineID = line[0:6].strip()                       # L     String 6
        if self.lineID == "ATOM" or self.lineID == "HETATM":
            self.atmNum = int(line[6:11])                     # R     Int 5
            self.atmName = line[12:16].strip()                # L     String 4
            if len(self.atmName) != 4:
                self.atmName = " " + self.atmName
            self.altLoc = line[16:17]                         #       String 1
            self.aaName = line[17:21].strip()                 # L     String 4
            self.chainID = line[21:22]                        #       String 1
            self.aaNum = int(line[22:26])                     # R     Int 4
            self.iCode = line[26:27]                          #       String 1
            self.x = float(line[30:38])                       # R     Float 8.3
            self.y = float(line[38:46])                       # R     Float 8.3
            self.z = float(line[46:54])                       # R     Float 8.3
	    try:					# If field is populated
            	self.occ = float(line[54:60])                 # R     Float 6.2
	    except:
		self.occ = 1.0
	    try:					# If field is populated
		self.tempFac = float(line[60:66])             # R     Float 6.2
	    except:
		self.tempFac = 0.0
            try:					# If field is populated
                self.ele = line[76:78]                        # R     String 2
            except:
                self.ele = ""
	    self.resatm = "%s.%s" %(self.aaName,self.atmName.strip())
        elif self.lineID == "TER":                      # Chain Termination
            self.info = line[6:].strip()
            self.atmName = line[12:16].strip()                # L     String 4
            self.altLoc = line[16:17]                         #       String 1
            self.aaName = line[17:21].strip()                 # L     String 4
            self.chainID = line[21:22]                        #       String 1
            try:                                        # If field is populated
                self.atmNum = int(line[6:11])                 # R     Int 5
                self.aaNum = int(line[22:26])                 # R     Int 4
            except:                                     # If field is populated
                self.atmNum = 0                               # R     Int 5
                self.aaNum = 0                                # R     Int 4
        elif self.lineID == "MODEL":                    # Model Number
            self.nModel = int(line[6:].strip())
        elif self.lineID == "DBREF":                    # Outside DB Reference
            self.ID = line[7:11] + ":" + line[12]       # pdbID:chainID
            self.chainID = line[12]                     # chainID
            self.unpID = line[33:41].strip()            # DatabaseID (often uniprot)
            self.nstart = int(line[14:18].strip())      # Start number of pdb sequence
            self.nend = int(line[20:24].strip())        # End number of pdb sequence
            self.start = line[55:60].strip()            # Start number of database sequence
            self.end = line[62:67].strip()              # End number of database sequence
        elif self.lineID == "HET":                      # Heteroatom record
            self.hetID = line[7:10].strip()             # HET molecule ID
            self.chainID = line[12]                     # HET molecule chainID
            self.hetNum = int(line[13:17].strip())      # HET molecule sequence number in chain
            self.iCode = line[17]
            self.numHetAtoms = int(line[20:25].strip()) # Number of associated HETATM records
        elif self.lineID == "MODRES":                   # Modified residue record
            self.hetID = line[12:15].strip()            # Modified residue name (HET or normal)
            self.chainID = line[16]                     # HET molecule chainID
            self.hetNum = int(line[18:22].strip())      # HET molecule sequence number in chain
	    self.iCode = line[22]
            self.aaName = line[24:27].strip()           # Standard residue name
        elif self.lineID == "SEQADV":                   # SEQADV record
            self.hetID = line[12:15].strip()            # PDB residue name (HET or normal)
            self.chainID = line[16]                     # PDB molecule chainID
            try:                                        # PDB molecule sequence number
                self.aaNum = int(line[18:22].strip())
            except:
                self.aaNum = line[18:22].strip()
            self.aaName = line[39:42].strip()           # Standard residue name
            self.note = line[49:].strip()               # Notes
            self.all = line
        elif self.lineID == "REMARK":                   # REMARKS
            if line[6:10].strip() == "610":             # REMARK 610
                self.aaName = line[15:18].strip()
                self.chainID = line[19]
                self.aaNum = int(line[21:25])
	    if line[6:10].strip() == "465":		# REMARK 465
		self.aaName = line[15:18].strip()
		self.chainID = line[19]
		self.aaNum = int(line[21:26])
		try:
		    self.iCode = line[26]
		except:
		    self.iCode = " "
	elif self.lineID == "SEQRES":
	    self.chainID = line[11]
	    self.aaNames = line[19:].split()
                                           
    def __str__(self):
        if self.lineID == "ATOM" or self.lineID == "HETATM":
            return ("%-6s%5d %-4s%1s%-4s%1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f"
                    "          %2s" %(self.lineID, self.atmNum, self.atmName, 
                    self.altLoc, self.aaName, self.chainID, self.aaNum, 
                    self.iCode, self.x, self.y, self.z, self.occ, self.tempFac,
                    self.ele))
        elif self.lineID == "TER":
            return ("%-6s%5d %-4s%1s%-4s%1s%4d" %(self.lineID, self.atmNum, 
                    self.atmName, self.altLoc, self.aaName, self.chainID, 
                    self.aaNum))
        elif self.lineID == "ENDMDL" or self.lineID == "END":
            return ("END")
	else:
	    return ("NO STR FUNCTION DEFINED")

    def distance(self,self2):
        dx = pow(self.x - self2.x, 2)
        dy = pow(self.y - self2.y, 2)
        dz = pow(self.z - self2.z, 2)
        return pow(dx + dy + dz, 0.5)
    
    def update_lineID(self,new_lineID):
        self.lineID = new_lineID

    def update_atmNum(self,new_atmNum):
        self.atmNum = new_atmNum

    def update_atmName(self,new_atmName):
        self.atmName = new_atmName
        if len(self.atmName) != 4:
            self.atmName = " " + self.atmName       # Pad a space in front
	self.resatm = "%s.%s" %(self.aaName,self.atmName.strip())

    def update_altLoc(self,new_altLoc):
        self.altLoc = new_altLoc
            
    def update_aaName(self,new_aaName):
        self.aaName = new_aaName
	self.resatm = "%s.%s" %(self.aaName,self.atmName.strip())

    def update_chainID(self,new_chainID):
        self.chainID = new_chainID

    def update_aaNum(self,new_aaNum):
        self.aaNum = new_aaNum

    def update_iCode(self,new_iCode):
        self.iCode = new_iCode
            
    def update_x(self,new_x):
        self.x = new_x

    def update_y(self,new_y):
        self.y = new_y

    def update_z(self,new_z):
        self.z = new_z

    def update_occ(self,new_occ):
        self.occ = new_occ

    def update_tempFac(self,new_tempFac):
        self.tempFac = new_tempFac

    def update_ele(self,new_ele):
        self.ele = new_ele
