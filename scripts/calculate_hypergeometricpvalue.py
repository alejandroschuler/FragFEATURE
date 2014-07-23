from argparse import ArgumentParser
from collections import Counter, defaultdict
from numpy import take
from scipy.stats import hypergeom
import itertools
import cPickle as pickle

from load_datafiles import load_resatm_KB_boundfrags
from supportingcode import group_ff_byresatm

###############################################################################
def main():
    parser = ArgumentParser(description="Calculate the hypergeometric p-value "
                            "for each fragment using fragment information from"
                            " k nearest neighbors (--kNN) and including only "
                            "fragments with sufficient information (--nFRAG)")
    parser.add_argument("--folder", type=str, dest="folder",
                        help="folder of input/output files")
    parser.add_argument("--kNN", type=int, nargs="+", dest="kNN",
                        help="number of nearest neighbors to use")
    parser.add_argument("--nFRAG", type=int, nargs = "+", dest="nFRAG",
                        help="minimum fragment information required")
    args = parser.parse_args()

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

    # Group the feature vectors by RES.ATM tag
    group_byresatm = group_ff_byresatm(result) 

    # Test all combinations of kNN and nFRAG
    params = itertools.product(args.kNN, args.nFRAG)
    for knn, nfrag in params:
        print ("Calculating fragment hypergeometric p-values <- "
               "kNN:%02d , nFrag:%03d" %(knn, nfrag))

        fragpvalues = defaultdict(lambda : defaultdict(dict))
        for resatm in group_byresatm:
            # Load the fragments bound by the resatm KB 
            fragdata = load_resatm_KB_boundfrags(resatm)

            for i in group_byresatm[resatm]:
                # Number of non-homologous neighbors (number of total balls)
		n_nonh = len(result[i]['NonH_ANNIDX'])

                # Collect the fragments bound by the k nearest neighbors
		observedfrags = []
		for j in range(0, knn):
		    idx = result[i]['NonH_ANNIDX'][j]
		    if fragdata[idx] != ['']:
			observedfrags.extend(fragdata[idx])

                # Retrieve the fragments bound by the non-homologous neighbors
		nonhfrags = take(fragdata,result[i]['NonH_ANNIDX'])
		nonhfrags = [item for sublist in nonhfrags for item in sublist]

                # Count the number of times each fragment is bound by 
                # non-homologous neighbors
		nonhfragcount = Counter(nonhfrags)

		for fragx in list(set(observedfrags)):
		    k = knn  			    # Number balls drawn
	    	    x = observedfrags.count(fragx)  # Number fragx balls drawn
		    m = nonhfragcount[fragx] 	    # Number fragx balls
		    # Calculate the hypergeometric pvalue (-log)
                    pvalue = -hypergeom.logsf(x-1, n_nonh, m, k)
                    fragpvalues[i][fragx] = pvalue

                # Remove rare fragments because not enough is known about them
		frags = fragpvalues[i].keys()
		for frag in frags:
		    if nonhfragcount[frag] < nfrag:
			del fragpvalues[i][frag]

        # Pickle the results
        # fragpvalues[i] = for entry i, dict of fragment pvalues (-log) for
        # the fragments bound by the k nearest non-homologous neighbors
        outfile = ("%s/nonhomologous.kNN%d.nFrag%d.hypergeo.pvar" 
                   %(args.folder, knn, nfrag))
        outfile = open(outfile, 'w')
        pickle.dump(dict(fragpvalues), outfile, protocol=2)
        outfile.close()


if __name__ == "__main__":
    main()
