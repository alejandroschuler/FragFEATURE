from argparse import ArgumentParser
from collections import defaultdict
import cPickle as pickle

from supportingcode import group_ff_byresatm
from load_datafiles import load_resatm_KB_ann
from load_datafiles import load_resatm_KB_boundfrags
###############################################################################
def main():
    parser = ArgumentParser(description="Prints the knowledge base members "
                            "whose fragent information was used to generate "
                            "the fragment prediction")
    parser.add_argument("--folder", type=str, dest="folder",
                        help="folder of input/output files")
    parser.add_argument("--kNN", type=int, dest="kNN",
                        help="number of nearest neighbors to use")
    parser.add_argument("--nFRAG", type=int, dest="nFRAG",
                        help="minimum fragment information required")
    parser.add_argument("--ffs", type=str, dest="ffs",
                        help="feature vectors invovled in the prediction")
    parser.add_argument("--fragID", type=int, dest="frag",
                        help="fragment of interest for this microenvironment set")
    args = parser.parse_args()

    args.ffs = map(int, args.ffs.split(","))

    # Load the nonhomologous nearest neighbors
    # result[i]['ANN']         = for entry i, annotation information
    # result[i]['NonH_dVEC']   = for entry i, list of dissimilarity to each
    # Residue.Atom KB member (non-homologous) (sorted by dissimilarity)
    # result[i]['NonH_ANNIDX'] = for entry i, list of annotation index to
    # Residue.Atom KB (non-homologous) (sorted by dissimilarity)
    infile = "%s/nonhomologous.pvar" %args.folder
    infile = open(infile, 'r')
    result = pickle.load(infile)
    infile.close()

    # Keep the feature vectors of interest
    keys = result.keys()
    for key in keys:
	if key not in args.ffs:
	    del result[key]

    # Group the feature vectors by residue types
    group_byresatm = group_ff_byresatm(result)

    print "_________Query_________\t",
    print "_____________Nearest Neighbors Information_____________"
    print "FF#\tFF-Label\tk\tPDBID\tFF-Label\tBound Ligand Atom(s)"

    # Track all pdbIDs associated with the nearest neighbors voting for the fragment
    all_pdbIDs = []
    # Track all ligands bound by the nearest neighbors voting for the fragment
    all_ligIDs = []

    # Process each feature vector of interest
    for resatm in group_byresatm:
	# Load the knowledge base
	anndata = load_resatm_KB_ann(resatm)
	fragdata = load_resatm_KB_boundfrags(resatm)

	for ff in group_byresatm[resatm]:
	    # Process each nearest neighbor index
	    for x in range(0, args.kNN):
	    	# Define the nearest neighbor index
	    	nnindex = result[ff]['NonH_ANNIDX'][x]

	    	# Get the meta data for this nearest neighbor
		record = anndata[nnindex]
		aainfo = record.aaInfo
	    	ligIDs = record.ligIDs
	    	ligtags = defaultdict(list)
		for entry in record.ligtags:
                    entry = entry.split(".")
                    ligtags[entry[0]].append(entry[1])

		# Check if predicted fragment is bound by the nearest neighbor
                # Mark the neighbor with an * if binding the fragment
	    	if args.frag in fragdata[nnindex]:
		    all_pdbIDs.append(record.tag)
		    all_ligIDs.extend(ligIDs)
		    flag = "*"
	    	else:
		    flag = " "

                text = "%s\t%s\t%d%s\t%s\t%s\t" %(ff, result[ff]['ANN'].info,
		                                  x, flag, record.tag, aainfo)
                for ligID in ligtags:
                    text += "%s:%s " %(ligID, ",".join(ligtags[ligID]))
                print text

    # Print out some summary information
    print
    print "Fragment:%d\t# FF-Vectors:%d" %(args.frag,len(args.ffs)) 
    print
    print "* = a nearest neighbor contributing to the fragment prediction"
    print "Summaries below only refer to the entries marked with *"
    print
    print "Summary: Ligand IDs"
    for ligID in set(all_ligIDs):
	print "%s:%d  " %(ligID,all_ligIDs.count(ligID)),
    print

    print "Summary: PDB IDs"
    for pdbID in set(all_pdbIDs):
	print "%s:%d  " %(pdbID,all_pdbIDs.count(pdbID)),
    print


if __name__ == "__main__":
    main()
