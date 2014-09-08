from subprocess import Popen,PIPE
from numpy import matrix, array
from linetypes import FF_point, PDB_line
import cPickle as pickle
from load_directorypaths import *

###############################################################################
# Given a command line command, call the command and return the output
def get_cmdline_output(cmd):
    data = Popen(cmd, shell=True, stdout=PIPE).communicate()[0]
    data = data.strip().split("\n")
    if data == ['']:
        return []
    else:
        return data


###############################################################################
# Read in PDB %ID clusters
# clustdict[pdbID.chainID] = cluster number
def load_seqclusterid(percentID):
    infile = "%s/clustIDs/clusters-%s.pvar" %(KB_HOME, percentID)
    if os.path.exists(infile):
	pass
    else:
	infile = "%s/clustIDs/bc-%s.pvar" %(KB_HOME, percentID)
    infile = open(infile, 'r')
    clustdict = pickle.load(infile)
    infile.close()
    return clustdict


###############################################################################
# Load the list of functional atoms
# fxndict[Residue.Atom] = distance cutoff
def load_standard_microenvironmentcenters():
    infile = "%s/fxn.atms.final.lst" %KB_HOME
    infile = open(infile, 'r')
    functionalatomdict = {}
    for line in infile:
        line = line.strip().split()
        functionalatom = line[0]
        dist = float(line[1])  # Distance cutoff b/t functional atom and ligand
        functionalatomdict[functionalatom] = dist
    infile.close()
    return functionalatomdict


###############################################################################
# Load the pseudo atom definitions
# fxncenters[Residue.Atom] = atoms comprising the pseudo atom location
def load_pseudoatom_microenvironmentcenters():
    infile = "%s/fxn.centers.lst" %KB_HOME
    infile = open(infile, 'r')
    fxncenters = {}
    for line in infile:
        line = line.strip().split("\t")
        res = line[0].split(".")[0]
        atms = line[1]
        fxncenters[res]=atms.split(" ")
    infile.close()
    return fxncenters


###############################################################################
# fxncutoffs[Residue.Atom] = max(distance between protein and ligand atom)
def load_fxnatm_cutoffs():
    infile = "%s/fxn.atms.cutoff.lst" %KB_HOME
    infile = open(infile, 'r')
    fxncutoffs = {}
    for line in infile:
	line = line.strip().split()
	resatom = line[0]
	cutoff = float(line[1])
	fxncutoffs[resatom] = cutoff
    return fxncutoffs

###############################################################################
# Load ff vector annotations for RES.ATM
def load_resatm_KB_ann(resatm):
    infile = "%s/%s/%s.ALL.annotation.txt" %(KB_HOME, resatm, resatm)
    infile = open(infile, 'r')
    anndataKB = map(FF_point, infile.readlines())
    infile.close()
    return anndataKB


###############################################################################
# Load ff vector property information for RES.ATM
def load_resatm_KB_prop(resatm, filter_flag=1):
    infile = "%s/%s/%s.ALL.property.pvar" %(KB_HOME, resatm, resatm)
    infile = open(infile, 'r')
    ffdataKB = pickle.load(infile)
    infile.close()
    if filter_flag:
        stdev = load_resatm_stdev(resatm)
        # Pull the indices of non-zero standard deviation features
        non0stdev_index = [index for index, element in enumerate(stdev) if element != 0]
        ffdataKB = ffdataKB[:,non0stdev_index]

    return ffdataKB


###############################################################################
# Load the fragments bound by each microenvironmnet of the resatm kB
def load_resatm_KB_boundfrags(resatm):
    infile = "%s/%s/%s.ALL.boundfrags.txt" %(KB_HOME, resatm, resatm)
    infile = open(infile, 'r')
    fragdata = infile.readlines()
    infile.close()
    for i in range(0, len(fragdata)):
        fragdata[i] = fragdata[i].strip().split("\t")
        if fragdata[i] != ['']:
            fragdata[i] = map(int, fragdata[i])
        else:
            fragdata[i] = []
    return fragdata


###############################################################################
# Load NR pre-calculated stdev to use for bit conversion for RES.ATM
def load_resatm_stdev(resatm):
    infile = "%s/%s/ID95/%s.ALL.nr.stdev.1.0" %(KB_HOME, resatm, resatm)
    infile = open(infile, 'r')
    stdev = infile.readlines()
    stdev = array(stdev).astype(float)
    infile.close()
    return stdev


###############################################################################
# Load ATOM atoms for a given pdbID
def load_atoms(filename):
    cmd = "grep -E '^ATOM' " + filename
    pdbatms = map(PDB_line, get_cmdline_output(cmd))
    return pdbatms


###############################################################################
# Load HETATM atoms for a given pdbID
def load_hetatoms(filename):
    cmd = "grep -E '^HETATM' " + filename
    hetatoms = map(PDB_line, get_cmdline_output(cmd))
    return hetatoms


###############################################################################
# Load ff vector annotation information
def load_ff_annotations(filename):
    infile = open(filename, 'r')
    anndata = map(FF_point, infile.readlines())
    infile.close()
    return anndata


###############################################################################
# Load ff vector property information
def load_ff_properties(filename):
    infile = open(filename, 'r')
    ffprop = pickle.load(infile)
    infile.close()
    return ffprop


###############################################################################
# Load the occurrence (count) of each fragment in a resatm KB
def load_resatm_fragment_counts(resatm):
    infile = "%s/%s/BoundFragCount.pvar" %(KB_HOME, resatm)
    infile = open(infile, 'r')
    fragcount = pickle.load(infile)
    infile.close()
    return fragcount












###############################################################################
# fragment_mapping[ligandID.ligandatom] = list of associated fragments
def load_fragment_mapping():
    infile = "%s/pdb_to_frag.pvar" %KB_HOME
    infile = open(infile, 'r')
    fragment_mapping = pickle.load(infile)
    infile.close()
    return fragment_mapping

###############################################################################
def load_frag_binding(frag_ID):
    fragdir = KB_HOME + '/f' + str(frag_ID)
    infile = fragdir + '/binding_list.pvar'
    infile = open(infile, 'r')
    micros_who_bind_querry = pickle.load(infile)
    infile.close()
    return micros_who_bind_querry
