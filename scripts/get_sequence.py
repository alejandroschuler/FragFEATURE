from argparse import ArgumentParser
from collections import defaultdict

from load_datafiles import load_atoms
from aminoacids import *

###############################################################################

def main():
    parser = ArgumentParser()
    parser.add_argument("--pdbID", type=str, dest="pdbID", 
                        help="PDB ID to analyze a pocket in")
    parser.add_argument("--folder", type=str, dest="folder",
                        help="folder of input/output files")
    args = parser.parse_args()
    ###########################################################################

    pdbID = args.pdbID.upper()

    # Load the pdb protein atoms
    pdbatms = load_atoms("%s/%s.pdb" %(args.folder, pdbID))
    pdbatms_bychain = defaultdict(list)
    for atm in pdbatms:
	if atm.aaName in threeLetterAA:  # Only include protien atoms
	    pdbatms_bychain[atm.chainID].append(atm)

    # Remove PDB chains that are peptides
    chainIDs = pdbatms_bychain.keys()
    for chainID in chainIDs:
        if len(pdbatms_bychain[chainID]) <= 20:
                del pdbatms_bychain[chainID]

    # Output the sequence of each protein chain
    for chainID in chainIDs:
	pdbatms = pdbatms_bychain[chainID]
        foundresidues = []  # Residues found
        for pdbatm in pdbatms:
            info = ["%s%s" %(pdbatm.aaNum,pdbatm.iCode),pdbatm.aaName]
            # If residue has not be seen, add it to foundresidues
            if info not in foundresidues:
                foundresidues.append(info)

        # Retrieve the amino acid sequence
        sequence = ""
        for entry in foundresidues:
            aa_3letter = entry[1]
            if aa_3letter in aa3to1dict:
                sequence += aa3to1dict[aa_3letter]
            else:
                sequence += "X"

	seqfile = "%s/%s.%s.fasta" %(args.folder, pdbID, chainID)
	seqfile = open(seqfile, 'w')
	seqfile.write(">%s.%s\n" %(pdbID, chainID))
	seqfile.write(sequence + "\n")
	seqfile.close()


if __name__ == "__main__":
    main()
