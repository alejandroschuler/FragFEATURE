import cPickle as pickle
import numpy  as np
import scipy.misc as sci
import random
from load_datafiles import load_resatm_KB_prop, load_resatm_stdev
from supportingcode import dissimilarity_to_KB
from load_directorypaths import *
from set_inclusion_core import *

def set_inclusion(querry_ID_string, test_resatms = [], test_num = 0):
### micros_who_bind_querry is a dictionary with residue atom types as keys 
### and lists of the line numbers of all that type of residue that are known 
### to bind the fragment as values

    random.seed(786914)

    infile = '%s/micros_who_bind_%s.pvar' % (KB_HOME, querry_ID_string)    
    infile = open(infile, 'r')
    micros_who_bind_querry = pickle.load(infile)

    # Determine the resatm types to test
    if test_resatms:
        print '\n\nChecking the requested microenvironment types'
        resatms = [resatm for resatm in test_resatms if resatm in micros_who_bind_querry.keys()]
        print 'Will calculate scores for: \n%s\n\n' % (str(resatms))
    else:
        resatms = micros_who_bind_querry.keys()

    # initialize the set inclusion score resatm dictionary
    intra_J_all_resatms = {}
    extra_J_all_resatms = {}

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
        J_list_resatm = compare_candidates_to_binding_set(
            intra_candidate_IDs, frag_binding_IDs, HomoFil, TanCalc)
        intra_J_all_resatms[resatm] = J_list_resatm

        J_list_resatm = compare_candidates_to_binding_set(
            extra_candidate_IDs, frag_binding_IDs, HomoFil, TanCalc)
        extra_J_all_resatms[resatm] = J_list_resatm

#    outfile = '%s/set_inclusion_scores_%s' %(KB_HOME, querry_ID_string)
#    outfile = open(outfile, 'w')
#    pickle.dump((J_intra,J_extra), outfile)

    return (intra_J_all_resatms, extra_J_all_resatms)
