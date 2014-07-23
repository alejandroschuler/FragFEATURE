from argparse import ArgumentParser
from collections import defaultdict
from numpy import unique, matrix, array, where, average
from linetypes import PDB_line
import sys

from load_datafiles import get_cmdline_output
from supportingcode import get_bound_ligands
from supportingcode import check_pdbfile_usability

###############################################################################
# Remove duplicate residues 
# Select altLoc with most number of atoms or highest average occupancy
def refineduplicates(data):
    if not data:
        return []

    keepatms = []  # Indices of atoms to keep
    data = matrix(data)

    # Compile unique list of residues
    reslabels = unique(array(data[:,1].flatten()))
    for label in reslabels:
        # Extract rows for this label
        x1,y1 = where(data == label)
        subdata1 = data[x1,:]
        # Compile unique list of altLoc flags
        altflags = unique(array(subdata1[:,-1].flatten()))

        altflags_count = {}  # Track the number of atoms per flag
        altflags_occ = {}    # Track the average occupancy of each flag
        for flag in altflags:
            # Extract rows for this label and flag
            x2,y2 = where(subdata1 == flag)
            altflags_count[len(x2)] = flag
            average_occ = average(subdata1[x2,2].astype(float).flatten())
            altflags_occ[average_occ] = flag

        if len(altflags_count.keys()) != 1:
            # Not all altLoc flags have same number of atoms (choose max atoms)
            keepflag = altflags_count[max(altflags_count.keys())]
        else:
            # All altLoc flags have same number of atoms (choose max occupancy)
            keepflag = altflags_occ[max(altflags_occ.keys())]

        x2,y2 = where(subdata1 == keepflag)
        keepatms.extend(subdata1[x2,0].astype(int).flatten().tolist()[0])
    return keepatms

###############################################################################
# Get atoms given a list of indices
def get_atoms(keepatms, pdbbody):
    data = []
    keepatms.sort()  # Sort atom indices in numerical order
    for idx in keepatms:
        line = pdbbody[idx]
        data.append(line)
    return data

###############################################################################
def pdbprint(folder, pdbdata, pdbbody_original):
    pdbID = pdbdata["pdbID"]
    pdbheader = pdbdata["pdbHeader"]
    hetheader = pdbdata["hetHeader"]

    # Remove duplicate atoms and merge with keep atoms
    pdbdata["keepAtms"].extend(refineduplicates(pdbdata["dupAtms"]))
    pdbdata["keepHetAtms"].extend(refineduplicates(pdbdata["dupHetAtms"]))

    # Get atom records using the indices
    pdbBody = get_atoms(pdbdata["keepAtms"], pdbbody_original)
    hetBody = get_atoms(pdbdata["keepHetAtms"], pdbbody_original)
    
    pdb = "%s/%s.pdb" %(folder, pdbID)
    pdb = open(pdb, 'w')
    pdb.write("%s\n%s\n" %('\n'.join(pdbheader), '\n'.join(map(str,pdbBody))))
    pdb.close()

    het = "%s/%s.het" %(folder, pdbID)
    het = open(het, 'w')
    het.write("%s\n%s\n" %('\n'.join(hetheader), '\n'.join(map(str,hetBody))))
    het.close()


###############################################################################
def main():
    parser = ArgumentParser()
    parser.add_argument("--file", type=str, dest="filename")
    parser.add_argument("--folder", type=str, dest="folder",
                        help="folder of input/output files")
    args = parser.parse_args()

    # Check if the PDB file is usable (no C-alpha models or UNK residues)
    if not check_pdbfile_usability(args.filename):
        return

    # Otherwise process the PDB file
    pdbdata = defaultdict(list)

    # Extract the PDB ID
    pdbdata["pdbID"] = args.filename.split("/")[-1][:-4]

    # Retrieve pdb file header information (list of strings)
    cmd  = "grep -E '^HEADER|^COMPND|^SOURCE|^AUTHOR|^EXPDTA|^SEQADV|"
    cmd += "^DBREF|^SEQRES|^MODRES|^HET   |^LINK|^SSBOND|^REMARK   2|"
    cmd += "^REMARK 465|^REMARK 470|^REMARK 480|^REMARK 610' "
    pdbdata["pdbHeader"] = get_cmdline_output(cmd + args.filename)

    # Retrieve hetatom header information (list of strings)
    cmd = "grep -E '^HET   |^HETNAM|^HETSYN|^FORMUL|^CONECT|^LINK' "
    pdbdata["hetHeader"] = get_cmdline_output(cmd + args.filename)

    # Retrieve PDB body (list of pdb_line objects)
    cmd = "grep -E '^MODEL[ ]+[0-9]|^ATOM|^HETATM|^TER|^END|^ENDMDL' "
    pdbbody = map(PDB_line, get_cmdline_output(cmd + args.filename))

    # Retrieve the ligand bound in the pdb file
    boundligands = get_bound_ligands(args.filename)

    i = 0   # Counter
    while (i < len(pdbbody)):
        line = pdbbody[i]
        # Atom record
        if line.lineID == "ATOM":
            # Keep CNOS (P for DNA) (X for ASX GLX)
            if line.ele.strip() in "CNOSPX":
	        # Keep records with no Alternate Location
                if line.altLoc == " ":
                    pdbdata["keepAtms"].append(i)
                    i = i + 1
		# Store duplicate atoms for processing
                else:
                    label = line.aaName + line.chainID + str(line.aaNum)
                    pdbdata["dupAtms"].append([i, label, line.occ, line.altLoc])
                    line.update_altLoc(" ")
                    i = i + 1
            # Exclude hydrogen records
            elif line.ele.strip() in "HD":
                i = i + 1
            else:
                print >> sys.stderr, "Error: Unexpected element " + line.ele
                i = i + 1

        # Hetero Atom record
        elif line.lineID == "HETATM":
            # Revert selenomethonine to methonine
            if line.aaName == "MSE":
                line.update_lineID("ATOM")
                line.update_aaName("MET")
                if line.atmName == " SE":
                    line.update_atmName("SD")
                if line.ele == "SE":
                    line.update_ele(" S")
            # Revert selenocysteine to cysteine
            elif line.aaName == "CSE":
                line.update_lineID("ATOM")
                line.update_aaName("CYS")
                if line.atmName == " SE":
                    line.update_atmName("SG")
                if line.ele == "SE":
                    line.update_ele(" S")
            else:
                hetlabel = "%s.%s.%s" %(line.aaName, line.aaNum, line.chainID)
                if hetlabel in boundligands:  # Accepted ligand
                    if line.ele.strip() in "HD":  # Remove hydrogens
                        i = i + 1
                    else:
                        # Keep records with no Alternate Location
                        if line.altLoc == " ":
                            pdbdata["keepHetAtms"].append(i)
                            i = i + 1
                        # Store duplicate atoms for processing
                        else:
                            label = line.aaName + line.chainID + str(line.aaNum)
                            pdbdata["dupHetAtms"].append([i, label, line.occ, 
                                                          line.altLoc])
                            line.update_altLoc(" ")
                            i = i + 1
                else:  # Not-accepted ligand
                    i = i + 1

	# Model start separator (keep only the first model)
	elif line.lineID == "MODEL":
	    if i == 0:
		i = i + 1
	    else:
		break

	# Chain separator or ATOM/HETATM separator
        elif line.lineID == "TER":
            pdbdata["keepAtms"].append(i)
            i = i + 1

	# Model end separator
        elif line.lineID == "ENDMDL" or line.lineID == "END":
            pdbdata["keepAtms"].append(i)
            i = i + 1

    # Output the cleaned pdb data
    pdbprint(args.folder, pdbdata, pdbbody)


if __name__ == "__main__":
    main()
