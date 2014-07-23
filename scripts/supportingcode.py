import cPickle as pickle
import numpy as np
import subprocess
import os
import copy
from collections import defaultdict

from load_directorypaths import *
from load_datafiles import get_cmdline_output, load_hetatoms
from load_datafiles import load_pseudoatom_microenvironmentcenters
from aminoacids import threeLetterAA
import het_excludes
import linetypes

###############################################################################
def define_resatm(annotation):
    # Define the resatms that are part of a larger group
    groupings = {}
    infile = "%s/fxn.atms.merge.lst" %KB_HOME
    infile = open(infile, 'r')
    for line in infile:
        line = line.strip().split()
        groupname = line[0]
        for resatm in line[1:]:
            groupings[resatm] = groupname
    infile.close()

    resatm = "%s.%s" %(annotation.aaName, annotation.atmName)
    if resatm in groupings:
	resatm = groupings[resatm]
    return resatm


###############################################################################
# Group the feature vectors by Residue.Atom
def group_ff_byresatm(data):
    group_byresatm = defaultdict(list)
    for i in data:
        resatm = define_resatm(data[i]['ANN'])
        group_byresatm[resatm].append(i)
    return group_byresatm


###############################################################################
# Make FF file from PTF
def make_ff_file(ptfname, ffname, logname, path):
    # Featurize the PTF file
    cmd = "featurize -P %s -s %s > %s 2> %s" %(ptfname, path, ffname, logname)
    subprocess.Popen(cmd, shell=True).wait()
    # Remove the comment lines from the FF file
    cmd = "grep -v ^\# %s > %s.temp" %(ffname, ffname)
    subprocess.Popen(cmd, shell=True).wait()
    os.rename(ffname + ".temp", ffname)


###############################################################################
# Process FF file into two parts, annotation and property values
def process_ff_file(ffname):
    fffile = open(ffname, 'r')
    ffdata = fffile.read().splitlines()
    fffile.close()

    ffann = [] ; ffprop = []
    for line in map(linetypes.FF_line, ffdata):
	ffann.append([line.ID, line.x, line.y, line.z, line.info])
	ffprop.append(line.props)
    ffprop = np.matrix(ffprop).astype(float)

    # Output FF vector property values
    ffproperties = ffname[:-3] + ".property.pvar"
    ffproperties = open(ffproperties, 'w')
    pickle.dump(ffprop, ffproperties, protocol=2)
    ffproperties.close()

    # Output FF vector annotation information
    ffannotation = ffname[:-3] + ".annotation.txt"
    ffannotation = open(ffannotation, 'w')
    for line in ffann:
       	ffannotation.write("\t".join(map(str, line)) + "\n")
    ffannotation.close()


#################################################################################
# Make PTF line from PDBATM
# hetAtms = ligand atoms in proximity to the pdbAtm
def make_PTFline(entry, pdbAtm, hetAtms, minDist):
    xyz_info = "%s\t%.3f\t%.3f\t%.3f\t#" %(entry, pdbAtm.x, pdbAtm.y, pdbAtm.z)
    res_info = "%s:%s:%s%s:%s" %(pdbAtm.aaName, pdbAtm.atmName.strip(), 
                                 pdbAtm.aaNum, pdbAtm.iCode.strip(), 
                                 pdbAtm.chainID)
    lig_info = "%d" %len(hetAtms)    # Number of proximal ligand atoms
    for hetAtm in hetAtms:           # Concatenate the proximal ligand atoms
        lig_info += "|%s:%s:%s:%s:%s" %(hetAtm.aaName, hetAtm.ele.strip(),
                                        hetAtm.atmName.strip(), hetAtm.aaNum,
                                        hetAtm.chainID)
    lig_info += "|%.3f" %(minDist)  # Note the min distance to the ligand atoms
    return "%s\t%s|%s\n" %(xyz_info, res_info, lig_info)


###############################################################################
# Given the cluster id of the query, make the nearest neighbors non-homologous
# to the query and to each other
#def make_nonhomologous(database_clusterids, clust_id, dissimilarityvector, indexvec):
def make_nonhomologous(database_clusterids, dissimilarityvector):
    # Get a unique set of clusterids for the database/knowledge base
    uniq_clusterids = list(set(database_clusterids))

    # clusterid = -1 reflects (a) a dummy value for no clusterID or 
    #                         (b) a protein that is homologous to the query
    if -1 in uniq_clusterids:
        del uniq_clusterids[uniq_clusterids.index(-1)]
    #if clust_id in uniq_clusterids:
    #    del uniq_clusterids[uniq_clusterids.index(clust_id)]

    # Find entries non-homologous to each other, keeping the closest entry
    database_clusterids = np.array(database_clusterids)
    nonh_dissimilarityvector = []
    nonh_annotationindex = []
    for clustid in uniq_clusterids:
        clustindex, = np.where(database_clusterids == clustid)
        indexmindist = clustindex[np.argmin(dissimilarityvector[clustindex])]
        nonh_dissimilarityvector.append(dissimilarityvector[indexmindist])
        #nonh_annotationindex.append(indexvec[indexmindist])
        nonh_annotationindex.append(indexmindist)

    # Sort the non-homologous KB entries distance
    sortorder = np.argsort(nonh_dissimilarityvector)
    nonh_dissimilarityvector = np.take(nonh_dissimilarityvector, sortorder)
    nonh_annotationindex = np.take(nonh_annotationindex, sortorder)
    return nonh_dissimilarityvector, nonh_annotationindex



###############################################################################
# Given the atoms whose center defines a pseudo atom, calculate the center
def calculate_pseudoatom(pdbatms):
    cx = []; cy = []; cz = []
    for pdbatm in pdbatms:
        cx.append(pdbatm.x)
        cy.append(pdbatm.y)
        cz.append(pdbatm.z)
    # Create a fake atom to store the information
    pseudoatm = copy.deepcopy(pdbatm)
    pseudoatm.update_atmName("PSEU")
    pseudoatm.update_atmNum(0)
    pseudoatm.update_ele("X")
    pseudoatm.update_x(round(np.average(cx), 3))
    pseudoatm.update_y(round(np.average(cy), 3))
    pseudoatm.update_z(round(np.average(cz), 3))
    return pseudoatm


###############################################################################
# Add the pseudo atoms given a set of atoms
def add_pseudoatoms(pdbatms):
    # Load the pseudo atom definitions
    # fxncenters[Residue.Atom] = atoms comprising the pseudo atom location
    fxncenters = load_pseudoatom_microenvironmentcenters()

    # Determine the residues requiring a pseudo atom
    # Group all atoms into unique residues
    group_byres = defaultdict(lambda : defaultdict(list))
    for atm in pdbatms:
        # For residues with a special functional center
        if atm.aaName in fxncenters:
            # For residue atoms used to define the special functional center
            if atm.atmName.strip() in fxncenters[atm.aaName]:
                reslabel = "%s.%s.%s" %(atm.aaNum, atm.chainID, atm.iCode)
                group_byres[atm.aaName][reslabel].append(atm)

    # Calculate the pseudo atoms
    for res in group_byres:
        for reslabel in group_byres[res]:
            # Check atoms needed to calculate the functional center are found
            if len(group_byres[res][reslabel]) == len(fxncenters[res]):
                pseudoatm = calculate_pseudoatom(group_byres[res][reslabel])
                pdbatms.append(pseudoatm)
    return


###############################################################################
# Calculate a microenvironment's dissimilarity to the knowledge base
# ffvec = feature vector with non-varying columns removed
# ffKB = matrix of feature vectors with non-varying columns removed
# non0stdev = property standard deviations with non-varying columns removed
def dissimilarity_to_KB(ffvec, ffKB, non0stdev):
    # Subtract the ffvec from the knowledge base
    diffmatrix = abs(ffKB - ffvec)
    # Convert the diffmatrix into 0s and 1s
    # Less than or equal to stdev (1); Greater than the stdev (0)
    bitmatrix = np.less(diffmatrix, non0stdev).astype(int)

    # Position with non-zero values in ffKB and ffvec (0/1 , 1/0 , 1/1)
    not00matrix = (ffKB != 0) + (ffvec != 0)

    # Numerator = number of positions not 0/0 within stdev of each other
    numer = np.sum(np.multiply(bitmatrix, not00matrix), axis=1)
    # Denominator = n Number of positions not 0/0
    denom = np.sum(not00matrix,axis=1)

    # Vector of tanimoto similarity to knowledge base
    tc_vec = np.rint(10000.0*numer/denom).astype('uint16')
    tc_vec = tc_vec.flatten().tolist()[0]
    # Convert similarity to dissimilarity
    d_vec  = 10000 - np.array(tc_vec)
    return d_vec


###############################################################################
# Determine if the pdbatm is within the distance cutoff of any of the hetatms
def proximal_to_hetatms(pdbatm, hetatms, cutoff):
    for hetatm in hetatms:
        d = round(pdbatm.distance(hetatm), 3)
        if d <= cutoff:
            return True
    return False

###############################################################################
# Determine the distance between the pdbatm and each of the hetatms
def distance_to_hetatms(pdbatm, hetatms):
    d_vec = []
    for hetatm in hetatms:
        d = round(pdbatm.distance(hetatm), 3)
        d_vec.append(d)
    return np.array(d_vec)


###############################################################################
# Check if a PDB file is usable 
# Not Usable: Model Type of CA ATOMS ONLY or UNKs (unknown amino acids) found
def check_pdbfile_usability(file):
    if get_cmdline_output("grep -E '^MDLTYP' %s | grep 'CA ATOMS ONLY'" %file):
	return False
    if get_cmdline_output("grep -E '^ATOM  ' %s | grep UNK" %file):
	return False
    return True

###############################################################################
# 
def get_bound_ligands(file):
    # Default list of bad ligands
    badhet = copy.deepcopy(het_excludes.excludes)

    # Find all bad ligands (Modified residues/Chromophores/Acetylation sites
    data = []
    data.extend(get_cmdline_output("grep -E '^MODRES' %s" %file))
    data.extend(get_cmdline_output("grep -E '^SEQADV' %s | grep -E "
                                   "'MODIFIED|CHROMOPHORE|ACETYLATION'" %file))
    if data:
        data = map(linetypes.PDB_line, data)
        for line in data:
            badhet.append(line.hetID)

    # Glycosylation sites have addtional ligand not directly attached
    # Sugars do not count in the presence of glycosylation sites
    data = get_cmdline_output("grep -E '^MODRES' %s | grep "
                              "'GLYCOSYLATION SITE'" %file)
    if data:
        badhet.extend(het_excludes.sugars)

    # Any het names found in the SEQRES are covalently linked to the protein
    data = get_cmdline_output("grep -E '^SEQRES' %s" %file)
    data = map(linetypes.PDB_line, data)
    threeLetterAA.extend(["MSE", "CSE"])
    for line in data:
        for aaName in line.aaNames:
            if aaName not in threeLetterAA:
                badhet.append(aaName)

    # Remove hetero atoms associated with bad ligands
    hetatms = load_hetatoms(file)
    deleteidx = []
    for i in range(0,len(hetatms)):
        if hetatms[i].aaName in badhet:
            deleteidx.append(i)
    hetatms = np.delete(hetatms,deleteidx)

    # Remove hetero atoms associated with incomplete molecules
    # Retrieve ligand instances with incomplete molecule information
    badhetatms = get_cmdline_output("grep -E '^REMARK 610' " + file)
    deleteidx = []
    if badhetatms:
        # First 6 lines are header information
        badhetatms = map(linetypes.PDB_line, badhetatms[6:])
        for badhetatm in badhetatms:
            badhetatm_tag = (badhetatm.aaName + badhetatm.chainID 
                            + str(badhetatm.aaNum))
            for i in range(0, len(hetatms)):
		hetatm_tag = (hetatms[i].aaName + hetatms[i].chainID 
                              + str(hetatms[i].aaNum))
                if badhetatm_tag == hetatm_tag:
                    deleteidx.append(i)
    hetatms = np.delete(hetatms,deleteidx)

    # Gather the ligand labels for the remaining "good" ligands
    hettags = []
    for i in range(0, len(hetatms)):
        hettags.append("%s.%s.%s" %(hetatms[i].aaName, hetatms[i].aaNum, 
                                    hetatms[i].chainID))
    hettags = np.unique(hettags)
    return hettags
