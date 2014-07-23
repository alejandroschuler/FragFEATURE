import os
from argparse import ArgumentParser
from numpy import take, where

from load_directorypaths import *
from load_datafiles import load_atoms, load_hetatoms
from load_datafiles import load_fxnatm_cutoffs
from supportingcode import proximal_to_hetatms, distance_to_hetatms
from supportingcode import make_PTFline

###############################################################################
def main():
    parser = ArgumentParser()
    parser.add_argument("--residue", type=str, dest="res",
                        help="target residue")
    parser.add_argument("--atom", type=str, dest="atom",
                        help="target functional atom")
    args = parser.parse_args()

    # Load the cutoff distance permitted between protein and ligands
    # fxncutoffs[Residue.Atom] = max(distance between protein and ligand atom)
    fxncutoffs = load_fxnatm_cutoffs()

    # Read list of pdb files to process
    pdblistfile = "%s/logs/pdb.list.res.prot.usable.het" %PDB_HOME
    pdblistfile = open(pdblistfile, 'r')
    pdblist = pdblistfile.read().splitlines()
    pdblistfile.close()

    # Setup a folder for the target functional atom
    fxnatm = "%s.%s" %(args.res, args.atom)
    datadir = "%s/%s" %(KB_HOME, fxnatm)
    if not os.path.exists(datadir):
	os.mkdir(datadir)

    # Setup PTF file for the target functional atom
    ptffile = "%s/%s.ALL.ptf" %(datadir, fxnatm)
    print ptffile
    ptffile = open(ptffile, 'w')

    # Find the functional atom in the pdb files
    for entry in pdblist:
	pdbID = entry.split()[0]
        folder = "%s/pdb/%s" %(PDB_HOME, pdbID[1:3])

	# Load pdb atoms proximal to the ligand
        pdbatms = load_atoms("%s/%s.proximal" %(folder, pdbID))

	# Load the hetero atoms
        hetatms = load_hetatoms("%s/%s.het" %(folder, pdbID))

	for pdbatm in pdbatms:
	    # Check pdbatm is a functional atom
	    if pdbatm.resatm == fxnatm:
		# Calculate the distance between pdbatm and all hetatms
		distvec = distance_to_hetatms(pdbatm, hetatms)
		# Distance threshold for this functional atom
		cutoff = fxncutoffs[fxnatm]

                # Retrieve all hetatms within the distance threshold
                proximalhets = take(hetatms, where(distvec <= cutoff))[0]
                text = make_PTFline(pdbID, pdbatm, proximalhets, min(distvec))
                ptffile.write(text)
    ptffile.close()


if __name__ == "__main__":
    main()
