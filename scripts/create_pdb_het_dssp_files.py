from argparse import ArgumentParser
from subprocess import Popen, PIPE
import os, sys

from supportingcode import check_pdbfile_usability

def main():
    ###########################################################################
    parser = ArgumentParser()
    parser.add_argument("--pdbID", type=str, dest="pdbID", 
                        help="PDB ID to analyze a pocket in")
    parser.add_argument("--folder", type=str, dest="folder",
                        help="folder of input/output files")
    args = parser.parse_args()
    ###########################################################################

    if not os.path.exists(args.folder):
	os.mkdir(args.folder)

    # Download the desired PDB file
    pdbID = args.pdbID.upper()
    cmd = "wget -q http://www.rcsb.org/pdb/files/%s.pdb" %pdbID
    process = Popen(cmd, shell=True)
    process.wait()

    if process.returncode != 0:  # PDB file failed to download
	print >> sys.stderr, "Download failed!"
	print >> sys.stderr, "Please download %s.pdb from the PDB" %pdbID
	return

    # Check if the PDB file is usable (no C-alpha models or UNK residues)
    if not check_pdbfile_usability("%s.pdb" %pdbID):
	print >> sys.stderr, ("Protein file is not usable (C-alpha model or "
                              "residues of unknown identity")
	return

    # Clean up the PDB file into a clean pdb file and hetatm file
    cmd = ("python scripts/clean_PDB.py --file %s.pdb --folder %s" 
           %(pdbID, args.folder))
    Popen(cmd, shell=True).wait()

    # Check if dssp program is installed
    process = Popen("which dsspcmbi", shell=True, stdout=PIPE, stderr=PIPE)
    stdout, stderr = process.communicate()
    if process.returncode == 0:  # program found
	# Create corresponding DSSP file for the pdb file
	cmd = ("dsspcmbi %s/%s.pdb %s/%s.dssp 2>> %s/%s.dsspcmbi.log" 
               %(args.folder, pdbID, args.folder, pdbID, args.folder, pdbID))
	Popen(cmd,shell=True).wait()
    else:  # program not found
        print >> sys.stderr, "Program:dsspcmbi not found in $PATH"
	print >> sys.stderr, ("Please install dsspcmbi or generate your own dssp file"
                              " for the protein")

    # Remove the PDB file downloaded
    os.remove("%s.pdb" %pdbID)


if __name__ == "__main__":
    main()
