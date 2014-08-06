import cPickle as pickle
import os
import random
from load_datafiles import load_resatm_KB_prop, load_resatm_stdev
from supportingcode import dissimilarity_to_KB
from load_directorypaths import *

from candidate_to_binding_tanimoto_core import *

def main(querry_ID, test_resatms = [], test_num = 0):
### micros_who_bind_querry is a dictionary with residue atom types as keys 
### and lists of the line numbers of all that type of residue that are known 
### to bind the fragment as values

    random.seed(786914)

    querry_ID_string = str(querry_ID)
    fragdir = KB_HOME + '/f' + querry_ID_string

    infile = fragdir + '/binding_list.pvar'  
    infile = open(infile, 'r')
    micros_who_bind_querry = pickle.load(infile)

    # Determine the resatm types to test
    if test_resatms:
        print '\n\nChecking the requested microenvironment types'
        resatms = [resatm for resatm in test_resatms if resatm in micros_who_bind_querry.keys()]
    else:
        resatms = micros_who_bind_querry.keys()
    resatms = sorted(resatms)
    print 'Will calculate scores for: \n%s\n\n' % (str(resatms))

    # initialize the Tanimoto score resatm dictionary
    T = dict.fromkeys(['intra','extra'])

    for resatm in resatms:
        
        # IDs of micros (of type resatm) that bind the querry fragment
        frag_binding_IDs = micros_who_bind_querry[resatm]
        # Objects for filtration of homologs and calculation of Tanimoto coefficients
        HomoFil, TanCalc = prep_set_data(resatm, frag_binding_IDs)
        # IDs of micros (of type resatm) that are not known to bind the querry fragment
        nonbinding_IDs = list( set(range(len(TanCalc.resatmKB))) - set(frag_binding_IDs) )

        # Determine the candidates (drawn from non-binding set) to compare to the frag binding set
        print 'DETERMINING CANDIDATE SETS\n'
        if test_num != 0:
            # Determine the number of candidates to test
            true_test_num = min([test_num, len(frag_binding_IDs), len(nonbinding_IDs)])
            # Grab a random subset of these micro IDs to be candidates
            intra_candidate_IDs = random.sample(frag_binding_IDs, true_test_num)
            extra_candidate_IDs = random.sample(nonbinding_IDs, true_test_num)
            print 'Will test the inclusion of %d candidates\n' % (true_test_num)
        else:
            intra_candidate_IDs = frag_binding_IDs
            extra_candidate_IDs = random.sample(nonbinding_IDs, len(frag_binding_IDs))
            print 'Will test the inclusion of all micros'
        
        # Calculate the inclusion scores of all candidates
        print 'CALCULATING INCLUSION SCORES FOR ALL CANDIDATES\n'
        T['intra'] = tanimoto_matrix(intra_candidate_IDs, frag_binding_IDs, HomoFil, TanCalc, resatm)
        T['extra'] = tanimoto_matrix(extra_candidate_IDs, frag_binding_IDs, HomoFil, TanCalc, resatm)

        outfile = '%s/%s_tanimotos_%dN.pvar' % (fragdir, resatm, true_test_num)
        outfile = open(outfile, 'w')
        pickle.dump(T, outfile)
        outfile.close()

if __name__ == "__main__":
    import sys
    main(sys.argv[1], )
