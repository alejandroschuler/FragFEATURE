from argparse import ArgumentParser

from load_datafiles import load_fxnatm_cutoffs
from load_datafiles import load_atoms, load_hetatoms
from supportingcode import add_pseudoatoms, proximal_to_hetatms

###############################################################################
def main():
    parser = ArgumentParser()
    parser.add_argument("--pdbID", type=str, dest="pdbID",
                        help="PDB ID to analyze a pocket in")
    parser.add_argument("--folder", type=str, dest="folder",
                        help="folder of input/output files")
    args = parser.parse_args()

    # Load the cutoff distance permitted between protein and ligands
    # fxncutoffs[Residue.Atom] = max(distance between protein and ligand atom)
    fxncutoffs = load_fxnatm_cutoffs()

    # Load the protein atoms and hetero atoms
    pdbatms = load_atoms("%s/%s.pdb" %(args.folder, args.pdbID))
    hetatms = load_hetatoms("%s/%s.het" %(args.folder, args.pdbID))

    # Calculate and add the pseudo atoms
    add_pseudoatoms(pdbatms)

    # Find the protein atoms in proximity to the ligands
    proximalatms = []
    for pdbatm in pdbatms:
        if pdbatm.resatm in fxncutoffs:    # Functional atom
            cutoff = fxncutoffs[pdbatm.resatm]
            # Keep if atom is in proximity to ligand atoms
            if proximal_to_hetatms(pdbatm, hetatms, cutoff):
                proximalatms.append(pdbatm)
            else: pass

    # Output the proximal atoms if they are found
    if proximalatms:
        outfile = "%s/%s.proximal" %(args.folder, args.pdbID)
        outfile = open(outfile, 'w')
        outfile.write("\n".join(map(str, proximalatms)) + "\n")
        outfile.close()


if __name__ == "__main__":
    main()
