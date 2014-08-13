import cPickle as pickle
import os
import random
from supportingcode import dissimilarity_to_KB
from load_directorypaths import *

from candidate_to_binding_tanimoto_core import *

def main(frag_ID, test_resatms = [], test_num = 0):
### micros_who_bind_querry is a dictionary with residue atom types as keys 
### and lists of the line numbers of all that type of residue that are known 
### to bind the fragment as values

   # random.seed(786914)

    fragdir = KB_HOME + '/f' + str(frag_ID)

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
            # Determine the number of candidates to test in each group
            true_test_num = min([test_num, len(frag_binding_IDs), len(nonbinding_IDs)])
            # Grab a random subset of these micro IDs to be candidates
            intra_candidate_IDs = random.sample(frag_binding_IDs, true_test_num)
            extra_candidate_IDs = random.sample(nonbinding_IDs, true_test_num)
            test_type = 'partial'
            print 'Will test the inclusion of %d candidates\n' % (true_test_num)
        else:
            true_test_num = min(2000, len(nonbinding_IDs))
            intra_candidate_IDs = frag_binding_IDs
            extra_candidate_IDs = random.sample(nonbinding_IDs, true_test_num)
            test_type = 'full'
            print 'Will test the inclusion of all micros'
        
        # Calculate the inclusion scores of all candidates
        print 'CALCULATING INCLUSION SCORES FOR ALL CANDIDATES\n'
        T_intra = tanimoto_matrix(intra_candidate_IDs, frag_binding_IDs, HomoFil, TanCalc, resatm)
        T_extra = tanimoto_matrix(extra_candidate_IDs, frag_binding_IDs, HomoFil, TanCalc, resatm)
        
        if T_intra.shape[1] != T_extra.shape[1]:
            print '!!!!!!!!! TANIMOTO MATRICIES NOT CORRECT DIMENSION !!!!!!'
            import sys
            sys.exit()

        outfile_intra = '%s/%s_intra_tanimotos_%s.pvar' % (fragdir, resatm, test_type)
        outfile_extra = '%s/%s_extra_tanimotos_%s.pvar' % (fragdir, resatm, test_type)

        T_intra.dump(outfile_intra)
        T_extra.dump(outfile_extra)

if __name__ == "__main__":
    import sys
    main(sys.argv[1])
