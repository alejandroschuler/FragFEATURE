from argparse import ArgumentParser
import cPickle as pickle

from load_datafiles import load_resatm_KB_prop, load_resatm_stdev
from supportingcode import group_ff_byresatm, dissimilarity_to_KB

###############################################################################
def main():
    parser = ArgumentParser(description="Calculate the dissimilarity between "
                            "each pocket microenvironment and the knowledge "
                            "base microenvironments")
    parser.add_argument("--folder", type=str, dest="folder",
                        help="folder of input/output files")
    parser.add_argument("--weight", type=float, dest="weight", default=1.0,
                        help="weight to apply to the property standard "
                        "deviation (default: 1.0)")
    args = parser.parse_args()

    # Load the microenvironments
    # result[i]['ANN']   = for entry i, annotation information
    # result[i]['FFVEC'] = for entry i, feature vector
    infile = "%s/microenvironments.pvar" %args.folder
    infile = open(infile, 'r')
    result = pickle.load(infile)
    infile.close()

    print "Comparing the pocket of interest to the knowledge base"
    # Group the feature vectors by ResidueAtom
    group_byresatm = group_ff_byresatm(result)

    for resatm in group_byresatm:
        # Load knowledge base for Residue.Atom
        resatmKB = load_resatm_KB_prop(resatm)

        # Load the pre-calculated standard deviation to use for bit conversion
        # for Residue.Atom
        stdev = load_resatm_stdev(resatm)

	# Keep just the non-zero variance properties
        resatmKB = resatmKB[:,(stdev != 0)]
        non0stdev = stdev[stdev != 0]

        for i in group_byresatm[resatm]:
            # Keep the feature vector properties with non-zero variance
            ffvec = result[i]['FFVEC'][:,(stdev != 0)]

            # Calculate dissimilarity to Residue.Atom knowledge base (KB)
            result[i]['dVEC'] = dissimilarity_to_KB(ffvec, resatmKB, non0stdev)
            del result[i]['FFVEC']

    # Pickle the results
    # result[i]['ANN']  = for entry i, annotation information
    # result[i]['dVEC'] = for entry i, list of distance to RES.ATM database
    outFile = "%s/dissimilaritymatrix.stdev%s.pvar" %(args.folder,args.weight)
    outFile = open(outFile,'w')
    pickle.dump(result,outFile,protocol=2)
    outFile.close()


if __name__ == "__main__":
    main()
