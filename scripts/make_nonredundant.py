from argparse import ArgumentParser
from collections import defaultdict
import os
import pickle

from load_directorypaths import *
from load_datafiles import load_seqclusterid, load_resatm_KB_prop
from load_datafiles import load_resatm_KB_ann

def main():
    ###########################################################################
    parser = ArgumentParser()
    parser.add_argument("--residue", type=str, dest="res", 
                        help="target residue")
    parser.add_argument("--atom", type=str, dest="atom", 
                        help="target functional atom")
    parser.add_argument("--nr-threshold", type=str, dest="id", 
                        help="sequence identity threshold for defining "
                             "'non-redundant'")
    args = parser.parse_args()
    ###########################################################################

    # Make a directory for the non-redundant feature vectors
    datadir = "%s/%s.%s/ID%s" %(KB_HOME, args.res, args.atom, args.id)
    if not os.path.exists(datadir):
	os.mkdir(datadir)

    # Read in resolution
    # resolution_dict[pdbID] = resolution
    infile = open("%s/resolu.pvar" %(KB_HOME), 'r')
    resolution_dict = pickle.load(infile)
    infile.close()

    # Read in PDB %ID clusters
    # clusterID_dict[pdbID.chainID] = cluster number
    clusterID_dict = load_seqclusterid(args.id)

    # Read in KB ff vector annotations
    anndata_KB = load_resatm_KB_ann("%s.%s" %(args.res, args.atom))

    # Group ff vectors by cluster ID and resolution
    group_byclusterID = defaultdict(lambda : defaultdict(list))
    for entry in anndata_KB:
	if entry.tag not in clusterID_dict: # ClusterID not assigned by the PDB
	    pass	
	else:
	    clusterID = clusterID_dict[entry.tag]
	    resolution = resolution_dict[entry.tag[:-2]]
	    group_byclusterID[clusterID][resolution].append(entry.tag)

    # For each cluster ID, select pdbID.chainID with best resolution
    nrlist = []	 # Gather th pdbID.chainIDs we want to keep
    for clusterID in group_byclusterID:
	minresolution = min(group_byclusterID[clusterID].keys())
	beststructure = group_byclusterID[clusterID][minresolution][0]
	nrlist.append(beststructure)

    # Determine the indices associated with the selected pdbID.chainID
    nrindices = []
    for i in range(0,len(anndata_KB)):
	if anndata_KB[i].tag in nrlist:
	    nrindices.append(i)

    # Output FF vector indices associated with non-redundant pdbID.chainID
    outfile = "%s/%s.%s.ALL.nr.idx.txt" %(datadir, args.res, args.atom)
    outfile = open(outfile, 'w')
    outfile.write("\n".join(map(str, nrindices)) + "\n")
    outfile.close()


if __name__ == "__main__":
    main()
