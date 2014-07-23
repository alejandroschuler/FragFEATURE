from argparse import ArgumentParser
from collections import defaultdict
from scipy.stats import chi2
import cPickle as pickle
import itertools

###############################################################################
def main():
    parser = ArgumentParser(description="")
    parser.add_argument("--folder", type=str, dest="folder",
                        help="folder of input/output files")
    parser.add_argument("--kNN", type=int, nargs="+", dest="kNN",
                        help="number of nearest neighbors to use")
    parser.add_argument("--nFRAG", type=int, nargs = "+", dest="nFRAG",
                        help="minimum fragment information required")
    args = parser.parse_args()

    # Load in the microenvironment sets
    infile = "%s/microenvironment.sets.txt" %args.folder
    infile = open(infile, 'r')
    setdata = infile.readlines()
    infile.close()

    ffsets = []
    for line in setdata:
        line = map(int, line.strip().split('\t')[0].split(','))
        ffsets.append(line)

    # Determine the kNN and nFRAG parameters to run
    params = itertools.product(args.kNN, args.nFRAG)

    # Run through the different parameter sets
    for knn,nfrag in params:
        print ("Calculating microenvironment set fragment p-values <- "
               "kNN:%02d , nFrag:%03d" %(knn, nfrag))

        # fragpvalues[i] = for entry i, dict of fragment pvalues (-log) for
        # the fragments bound by the k nearest non-homologous neighbors
        infile = ("%s/nonhomologous.kNN%d.nFrag%d.hypergeo.pvar" 
                  %(args.folder, knn, nfrag))
        infile = open(infile, 'r')
        fragpvalues = pickle.load(infile)
        infile.close()

        # Store the results of aggregated feature vector predictions
        aggregate = defaultdict(dict)

        # Group feature vector predictions by their assignment
        for r in range(0, len(ffsets)):
            N = len(ffsets[r])

            # Aggregate fragment pvalues and occurrence for each set of 
            # microenvironments
            # Initialize the probability (-log) of observing a fragment to 0
            # Initialize the number of microenvironments predicting a fragment to 0
            set_fragpvalues = defaultdict(lambda: 0.0)
            set_fragCount = defaultdict(lambda: 0)
            for idx in ffsets[r]:
                for frag in fragpvalues[idx]:
                    set_fragpvalues[frag] += fragpvalues[idx][frag]
                    set_fragCount[frag] += 1

            # Delete fragments that are not voted for by all members of the set
            frags = set_fragpvalues.keys()
            for frag in frags:
                if set_fragCount[frag] != N:
                    del set_fragpvalues[frag]
                    del set_fragCount[frag]

            # Convert fragment pvalues to Fisher's chi-squared pvalue
            for frag in set_fragpvalues:
                set_fragpvalues[frag] = chi2.sf(2*set_fragpvalues[frag],2*N)

            if set_fragpvalues:  # Fragments predicted (not empty)
                aggregate[r] = dict(set_fragpvalues)

        # Output the remaining sets
        # aggregate[r] = for set r, dict of fragment pvalues 
        outfile = ("%s/nonhomologous.kNN%d.nFrag%d.fishers.pvar" 
                   %(args.folder, knn, nfrag))
        outfile = open(outfile, 'w')
        pickle.dump(aggregate, outfile, protocol=2)
        outfile.close()


if __name__ == "__main__":
    main()
