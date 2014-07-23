from argparse import ArgumentParser
from subprocess import Popen, PIPE
from collections import defaultdict
import os, sys, shutil

from load_datafiles import load_fxnatm_cutoffs
from load_datafiles import load_hetatoms, load_atoms, get_cmdline_output
from supportingcode import proximal_to_hetatms, add_pseudoatoms
from aminoacids import threeLetterAA

###############################################################################
def pocket_boundligand(folder, pdbID, pdbatms_bychain, fxncutoffs, result):
    # Load the hetero atoms for this pdb file & group by ligand instance
    hetatms = load_hetatoms("%s/%s.het" %(folder, pdbID))
    groupHetByLigTag = defaultdict(list)
    for atm in hetatms:
        ligTag = "%s.%s.%s" %(atm.aaName, atm.aaNum, atm.chainID)
        groupHetByLigTag[ligTag].append(atm)
    ligandTags = groupHetByLigTag.keys()
    ligandTags.sort()

    # Find functional atoms involved in binding each ligand instance
    for ligandTag in ligandTags:
        datadir = "%s/pocket_%s" %(folder, ligandTag)

        hetatms = groupHetByLigTag[ligandTag]  # Ligand atoms
        chainID = ligandTag.split(".")[-1]     # Ligand chainID
        pdbatms = pdbatms_bychain[chainID]     # pdbatms for this chainID

        for pdbatm in pdbatms:
            if pdbatm.resatm in fxncutoffs:    # Functional atom
                cutoff = fxncutoffs[pdbatm.resatm]
                if proximal_to_hetatms(pdbatm, hetatms, cutoff):
                    result[datadir].append(pdbatm)
                else: pass

def pocket_fpocket(folder, pdbID, pdbatms_bychain, fxncutoffs, maxpockets, result):
    for chainID in pdbatms_bychain:
	#  Write out each PDB protein chain (pseudo atoms do not count)
        pdbatms = pdbatms_bychain[chainID]
        outName = "%s/%s.%s.pdb" %(folder, pdbID, chainID)
        outFile = open(outName,'w')
        for pdbAtm in pdbatms:
	    if pdbAtm.atmName != "PSEU":
                outFile.write(str(pdbAtm) + "\n")
        outFile.close()

	# Find the pockets within the chain
        cmd = "fpocket -f %s" %(outName)
        pocketDir = "%s/%s.%s_out" %(folder, pdbID, chainID)
        Popen(cmd,shell=True).wait()
	os.remove(outName)

	numpockets = 0
        if os.path.exists(pocketDir):      # Pockets found
            cmd = "grep Pocket %s/%s.%s_info.txt | wc -l" %(pocketDir, pdbID, chainID)
            numpockets = int(get_cmdline_output(cmd)[0])
        numpockets = min([numpockets, maxpockets])

        datadir = "%s/pocket_f%s.%s" %(folder, maxpockets, chainID)

	# Retrieve the protein atoms associated with pockets (up to numpockets)
	pocketAtms = []
        for p in range(0,numpockets):
            pocket = "%s/pockets/pocket%d_atm.pdb" %(pocketDir,p)
            pocketAtms.extend(load_atoms(pocket))
        shutil.rmtree(pocketDir)

        # Group the pocket atoms into unique residues
        groupByResTag = defaultdict(list)
        for pdbAtm in pocketAtms:
            resTag = "%s.%s.%s" %(pdbAtm.aaName,pdbAtm.aaNum,pdbAtm.iCode)
            groupByResTag[resTag].append(pdbAtm)

        # For each residue, check if the pocket atoms correspond to residue sidechains
        sidechainResTags = []
        for resTag in groupByResTag:
            for pdbAtm in groupByResTag[resTag]:
		# Proline backbone is part of sidechain
                if pdbAtm.aaName == "PRO" and pdbAtm.atmName.strip() in ["N","CA"]:
                    sidechainResTags.append(resTag)
		# All other backbone atoms are backbones
                elif pdbAtm.atmName.strip() in ["N","C","O","CA"]:
                    if pdbAtm.resatm in fxncutoffs:  # Functional Backbone Atom
                        result[datadir].append(pdbAtm)
                else:
                    sidechainResTags.append(resTag)

        # Retrieve the functional atoms for the sidechains
        for pdbAtm in pdbatms:
            resTag = "%s.%s.%s" %(pdbAtm.aaName,pdbAtm.aaNum,pdbAtm.iCode)
            if resTag in sidechainResTags:
                if (pdbAtm.resatm in fxncutoffs) and (pdbAtm.atmName.strip() not in ["N","C","O","CA"]):
                    result[datadir].append(pdbAtm)


def pocket_residuelist(folder, pdbatms_bychain, fxncutoffs, residuefile, result):
    # Retrieve list of residue numbers for the pocket
    residuefile = "%s/%s" %(folder, residuefile)
    residuefile = open(residuefile, 'r')
    residuelist = defaultdict(list)
    for line in residuefile:
	line = line.split()
	resnum = int(line[0])
	chainID = line[1]
	residuelist[chainID].append(resnum)
    residuefile.close()

    # Find the protein atoms associated with the pocket
    for chainID in pdbatms_bychain:
        pdbatms = pdbatms_bychain[chainID]

	datadir = "%s/pocket_reslist.%s" %(folder, chainID)

        for atm in pdbatms:
            if atm.aaNum in residuelist[chainID]:
                if atm.resatm in fxncutoffs:   # Functional atom
                    result[datadir].append(atm)

###############################################################################

def main():
    parser = ArgumentParser()
    parser.add_argument("--pdbID", type=str, dest="pdbID", 
                        help="PDB ID to analyze a pocket in")
    parser.add_argument("--folder", type=str, dest="folder",
                        help="folder of input/output files")
    parser.add_argument("--boundligand", action="store_const", default=0, 
                        const=1, dest="pocket_mode", help="Use the bound "
                        "ligands to identify the pockets of interest")

    parser.add_argument("--fpocket", action="store_const", default=0, const=2,
                        dest="pocket_mode", help="Use fpocket to identify the"
                        " pockets of interest; must be used with --maxpockets")
    parser.add_argument("--maxpockets", type=int, dest="maxpockets",
                        help="Maximum number of pockets to use from fpocket; "
                        "must be used with --fpocket")

    parser.add_argument("--residuelist", action="store_const", default=0, 
                        const=3, dest="pocket_mode", help="Use a list of "
                        "residues to identify the pocket of interest; must be "
                        "used with --residuefile")
    parser.add_argument("--residuefile", type=str, dest="residuefile",
                        help="file with the list of residues comprising the "
                             "protein pocket; must be used with --fpocket")
    args = parser.parse_args()
    ###########################################################################

    # Load the cutoff distance permitted between protein and ligands
    # fxncutoffs[Residue.Atom] = max(distance between protein and ligand atom)
    fxncutoffs = load_fxnatm_cutoffs()

    # Load the pdb protein atoms
    pdbatms = load_atoms("%s/%s.pdb" %(args.folder, args.pdbID))
    pdbatms_bychain = defaultdict(list)
    for atm in pdbatms:
	if atm.aaName in threeLetterAA:  # Only include protien atoms
	    pdbatms_bychain[atm.chainID].append(atm)

    # Remove PDB chains that are peptides
    chainIDs = pdbatms_bychain.keys()
    for chainID in chainIDs:
        if len(pdbatms_bychain[chainID]) <= 20:
                del pdbatms_bychain[chainID]

    # Calculate the pseudoatoms
    chainIDs = pdbatms_bychain.keys()
    for chainID in chainIDs:
        pdbatms = pdbatms_bychain[chainID]
        add_pseudoatoms(pdbatms)

    proximalatms = defaultdict(list)
    if args.pocket_mode in [0, 1]:
   	pocket_boundligand(args.folder, args.pdbID, pdbatms_bychain, fxncutoffs, proximalatms)
    if args.pocket_mode in [0, 2]:
	process = Popen("which fpocket", shell=True, stdout=PIPE, stderr=PIPE)
	stdout, stderr = process.communicate()
	if process.returncode != 0:  # program not found
	    print >> sys.stderr, "Program:fpocket not found in $PATH"
	elif not args.maxpockets:
	    print >> sys.stderr, "Please specify a value for --maxpockets"
	else:
    	    pocket_fpocket(args.folder, args.pdbID, pdbatms_bychain, fxncutoffs, args.maxpockets, proximalatms)
    if args.pocket_mode in [0, 3]:
	if not args.residuefile:
	    print >> sys.stderr, "Please specify a value for --residuefile"
	else:
	    pocket_residuelist(args.folder, pdbatms_bychain, fxncutoffs, args.residuefile, proximalatms)

    # Output proximal atoms for the ligand
    for datadir in proximalatms:
        if not os.path.exists(datadir):
            os.mkdir(datadir)

	# Retrieve the labels for these proximal atoms
	output = []	
	for atm in proximalatms[datadir]:
	    output.append("%s:%s:%s:%s\n" %(atm.aaName, atm.atmName.strip(),
                                            atm.aaNum, atm.chainID))
	# Make sure duplicate atoms are only printed once
	# Duplicates should only be observed with args.pocket_mode = 2
	# Pockets can be overlapping, leading to duplicate atoms
	output = list(set(output))
	output.sort()

	# Output the atoms comprising this protein pocket
        outfile = "%s/%s.atoms" %(datadir, args.pdbID)
	outfile = open(outfile, 'w')
	outfile.write("".join(output))
	outfile.close()

	# Copy the pdb file and dssp file into the pocket folder
	shutil.copy2("%s/%s.pdb" %(args.folder, args.pdbID), datadir)
	shutil.copy2("%s/%s.dssp" %(args.folder, args.pdbID), datadir)

    pockets = proximalatms.keys()
    pockets.sort()
    print "Protein Pockets"
    print "\n".join(pockets)


if __name__ == "__main__":
    main()
