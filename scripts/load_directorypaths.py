import os
###############################################################################
# Define the necessary directories relating to running FragFEATURE
# Location of FragFEATURE
FragFEATURE_HOME = os.environ["HOME"] + "/FragFEATURE"

# Location of FragFEATURE Knowledge Base
KB_HOME	= FragFEATURE_HOME + "/knowledgeBase"

###############################################################################
# Define the necessary directories relating to creatng the Knowledge Base
PDB_HOME = FragFEATURE_HOME + "/sourcefiles"
# Location of pdb files listing atoms within 5 angstroms of a ligand
#PDB_HOME = "~/dbPDB"
#PDB_HET = PDB_HOME + "/het"
#PDB_WIN = PDB_HOME + "/within5"
