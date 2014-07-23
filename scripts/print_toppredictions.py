from argparse import ArgumentParser
from collections import defaultdict
import cPickle as pickle
import itertools
import os

###############################################################################
def main():
    parser = ArgumentParser(description="Prints the most significant fragment "
                            "predicted for a microenvironment set that meets "
                            "the significance threshold defined by --pvalue")
    parser.add_argument("--folder", type=str, dest="folder",
                        help="folder of input/output files")
    parser.add_argument("--kNN", type=int, nargs="+", dest="kNN",
                        help="number of nearest neighbors to use")
    parser.add_argument("--nFRAG", type=int, nargs = "+", dest="nFRAG",
                        help="minimum fragment information required")
    parser.add_argument("--pvalue", type=float, dest="pvalcutoff",
                        help="threshold p-value for significance")
    args = parser.parse_args()

    # Setup a directory for performance output
    if not os.path.exists("%s/predictions" %args.folder):
        os.mkdir("%s/predictions" %args.folder)

    # Determine the kNN and nFRAG parameters to run
    params = itertools.product(args.kNN, args.nFRAG)

    for knn,nfrag in params:
        print "Printing predictions <- kNN:%02d , nFrag:%03d" %(knn, nfrag)

        # Load the fragment preferences of the microenvironment sets
	# result[r] = for set r, dict of fragment pvalues
        infile = ("%s/nonhomologous.kNN%d.nFrag%d.fishers.pvar" 
                  %(args.folder, knn, nfrag))
        infile = open(infile, 'r')
        result = pickle.load(infile)
        infile.close()

        # Load in the microenvironment sets
        infile = "%s/microenvironment.sets.txt" %args.folder
        infile = open(infile, 'r')
        setdata = infile.readlines()
        infile.close()

        details = defaultdict(lambda : defaultdict(list))
        for i in range(0,len(setdata)):
            line = setdata[i].split()
            details[i]['FFIDX'] = map(int, line[0].split(','))
            details[i]['SETSPREAD'] = float(line[1])
            details[i]['ATOMS'] = line[2]
            details[i]['NN-SIMILARITY'] = line[3]

	# Remove the non-significant sets and predictions
        setIDs = result.keys()
        for idx in setIDs:
            pvals = result[idx].values()
            if min(pvals) >= args.pvalcutoff:
                del result[idx]
            else:
                frags = result[idx].keys()
                for frag in frags:
                    if result[idx][frag] >= args.pvalcutoff:
                        del result[idx][frag]

        # Absorb smaller sets into larger sets if smaller set is not as significant
        setIDs = result.keys()
        setIDs.sort()
        for i in range(0, len(setIDs)):
            set1idx = set(details[setIDs[i]]['FFIDX'])
            set1pval = min(result[setIDs[i]].values())
            for j in range(i+1,len(setIDs)):
                set2idx = set(details[setIDs[j]]['FFIDX'])
                set2pval = min(result[setIDs[j]].values())
                if set1idx.issubset(set2idx) and set1pval > set2pval:
                    del result[setIDs[i]]
                    break

        # Output the prediction results
        outfile = ("%s/predictions/nonhomologous.NN%d.nFrag%d.%.0e.predictions"
                   ".txt" %(args.folder, knn, 
                    nfrag, args.pvalcutoff))
        outfile = open(outfile, 'w')
        outfile.write("Fragment\tP-Value \tFF-IDs  \tSpread"
                      "\tNearest-Neighbor-Dissimilarity\tFunctional-Atoms\n")
        for i in result.keys():
            setfrag = min(result[i], key=result[i].get)
            setpval = min(result[i].values())

            outfile.write("%-10s\t%.2e\t" %(setfrag, setpval))
            outfile.write("%-8s\t" %(",".join(map(str, details[i]['FFIDX']))))
            outfile.write("%.2f\t" %(details[i]['SETSPREAD']))
            outfile.write("%-25s\t" %(details[i]['NN-SIMILARITY']))
            outfile.write("%s\n" %(details[i]['ATOMS']))

        outfile.close()


if __name__ == "__main__":
    main()
