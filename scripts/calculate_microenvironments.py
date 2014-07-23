from argparse import ArgumentParser
from collections import defaultdict
import cPickle as pickle
import sys
import os

from load_datafiles import load_pseudoatom_microenvironmentcenters, load_atoms
from load_datafiles import load_ff_annotations, load_ff_properties
from load_datafiles import get_cmdline_output
from supportingcode import make_ff_file, process_ff_file, calculate_pseudoatom
###############################################################################
def main():
    parser = ArgumentParser(description="Calculate microenvironments given a "
                            "list of atoms")
    parser.add_argument("--folder", type=str, dest="folder",
                        help="folder of input/output files")
    parser.add_argument("--file", type=str, dest="file",
                        help="file containing the atoms that define the pocket"
                             " of interest")
    args = parser.parse_args()

    # Define the name of the protein
    filename = ".".join(args.file.split(".")[:-1])

    # Load the pseudo atom definitions
    # fxncenters[Residue.Atom] = atoms comprising the pseudo atom location
    fxncenters = load_pseudoatom_microenvironmentcenters()

    # Load the atoms comprising the pocket of interest
    infile = "%s/%s" %(args.folder, args.file)
    infile = open(infile, 'r')
    pocketatoms = infile.read().splitlines()
    infile.close()

    # Load the atoms of the PDB file
    pdbatms = load_atoms("%s/%s.pdb" %(args.folder, filename))

    # Determine the residues requiring a pseudo atom by grouping all atoms for 
    # a residue into a dictionary entry
    group_byresidue = defaultdict(lambda : defaultdict(list))				
    for pdbatm in pdbatms:
        # If a residue has a special functional center, keep the atoms 
        # associated with the functional center definition
        if pdbatm.aaName in fxncenters:
            if pdbatm.atmName.strip() in fxncenters[pdbatm.aaName]:
                residuetag = "%s.%s.%s" %(pdbatm.aaNum, pdbatm.chainID, 
                                          pdbatm.iCode)
                group_byresidue[pdbatm.aaName][residuetag].append(pdbatm)

    # Calculate the x,y,z coordinate of the pseudo atoms and append these to 
    # the pdbatms variable
    for res in group_byresidue:
        for residuetag in group_byresidue[res]:
            # Ensure all atoms needed to calculate a functional center are present
            if len(group_byresidue[res][residuetag]) == len(fxncenters[res]):
                pseudoAtm = calculate_pseudoatom(group_byresidue[res][residuetag])
                pdbatms.append(pseudoAtm)

    # Output the pdb atoms that are part of the pocket in PTF file format
    ptf_filename = "%s/%s.ptf" %(args.folder, filename)
    outfile = open(ptf_filename, 'w')
    for pdbatm in pdbatms:
        atomtag = "%s:%s:%s:%s" %(pdbatm.aaName, pdbatm.atmName.strip(),
                                  pdbatm.aaNum, pdbatm.chainID)
        if atomtag in pocketatoms:
            outfile.write("%s\t%.3f\t%.3f\t%.3f\t#\t%s:%s:%s:%s\n" %(filename,
                          pdbatm.x, pdbatm.y, pdbatm.z, pdbatm.aaName,
                          pdbatm.atmName.strip(), pdbatm.aaNum, pdbatm.chainID))
    outfile.close()

    # Make FF file
    ff_logfile = "%s/%s.featurize.log" %(args.folder, filename)
    ff_filename = "%s/%s.ff" %(args.folder, filename)
    make_ff_file(ptf_filename, ff_filename, ff_logfile, args.folder)

    # Process FF file into two parts, annotation and property values
    process_ff_file(ff_filename)

    # Load FF annotation and property information
    annfile = ff_filename[:-3] + ".annotation.txt"
    anndata = load_ff_annotations(annfile)
    propfile = ff_filename[:-3] + ".property.pvar"
    propdata = load_ff_properties(propfile)

    # Repackage into a dictionary
    result = {}
    count = 0
    for i in range(0, len(anndata)):
        result[count] = {}
        result[count]['ANN'] = anndata[i]
        result[count]['FFVEC'] = propdata[i]
        count = count + 1

    # Pickle the results
    # result[i]['ANN']   = for entry i, annotation information
    # result[i]['FFVEC'] = for entry i, feature vector
    outfile = "%s/microenvironments.pvar" %args.folder
    outfile = open(outfile,'w')
    pickle.dump(result,outfile,protocol=2)
    outfile.close()

    # Check the log file for errors
    errors = get_cmdline_output("grep Error %s" %ff_logfile)
    if errors:
        print >> sys.stderr, errors
    else:
        os.remove(ff_logfile)

    # Remove files
    os.remove(ptf_filename)
    os.remove(ff_filename)
    os.remove(annfile)
    os.remove(propfile)

    print "Pocket of interest has %d microenvironments" %len(result)


if __name__ == "__main__":
    main()
