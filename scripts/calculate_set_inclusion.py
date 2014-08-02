import cPickle as pickle
import numpy  as np
import scipy.misc as sci
import random
from load_datafiles import load_resatm_KB_prop, load_resatm_stdev
from supportingcode import dissimilarity_to_KB
from load_directorypaths import *
from set_inclusion_core import *

def main(querry_ID_string, test_resatms=[], test_num=0):

    infile = '%s/micros_who_bind_%s.pvar' % (KB_HOME, querry_ID_string)    
    infile = open(infile, 'r')
    micros_who_bind_querry = pickle.load(infile)

    J_intra = intra_set_sim(micros_who_bind_querry, test_resatms, test_num)
    J_extra = extra_set_sim(micros_who_bind_querry, test_resatms, test_num)

    outfile = '%s/intra_set_inclusion_scores_%s' %(KB_HOME, querry_ID_string)
    outfile = open(outfile, 'w')
    pickle.dump(J, outfile)
    
    return (J_intra, J_extra)

#############################################################################3#

def intra_set_sim(micros_who_bind_querry, test_resatms = [], test_num = 0):
### micros_who_bind_querry is a dictionary with residue atom types as keys 
### and lists of the line numbers of all that type of residue that are known 
### to bind the fragment as values

    print 'CALCULATING SET INCLUSION FOR FRAG-BINDING MICRO-ENVIRONMENTS (EXPERIMENT)'

    if test_resatms:
        print 'Checking the requested microenvironment types'
        resatms = [resatm for resatm in test_resatms if resatm in micros_who_bind_querry.keys()]
        print 'Will calculate scores for: \n%s' % (str(resatms))
    else:
        resatms = micros_who_bind_querry.keys()

    # initialize the set inclusion score resatm dictionary
    J_all_resatms = {}

    for resatm in resatms:
        
        # Prepare the data for this resatm type
        print 'PRE-PROCESSING DATA FOR CALCULATION\n'
        frag_binding_ff_dict, non0stdev, resatm_clusters, resatmKB = (
            prep_set_data(resatm, micros_who_bind_querry[resatm]))

        print 'DETERMINING CANDIDATE SET\n'
        # Determine the candidates (drawn from non-binding set) to compare to the frag binding set
        if test_num != 0:
            # Determine the number of candidates to test
            true_test_num = min(test_num, len(frag_binding_ff_dict))
            # Grab a random subset of these micro IDs to be candidates
            candidate_ID_set = random.sample(frag_binding_ff_dict.keys(), true_test_num)
            # Generate a dictionary of ID# : ff_vector for the candidate microenvironments
            candidate_ff_dict = {candidate_ID : frag_binding_ff_dict[candidate_ID] for candidate_ID in candidate_ID_set}
            print 'Will test the inclusion of %d candidates\n' % (len(candidate_ff_dict))
        else:
            candidate_ff_dict = {candidate_ID : resatmKB[candidate_ID] for candidate_ID in nonbinding_ID_set} 
            print 'Will test the inclusion of entire frag-binding set (%d candidates)' % (len(candidate_ff_dict))
        
        # Calculate the inclusion scores of all candidates
        print 'CALCULATING INCLUSION SCORES FOR ALL CANDIDATES\n'
        J_list_resatm = (
            compare_candidates_to_binding_set(candidate_ff_dict, frag_binding_ff_dict, non0stdev, resatm_clusters))
        
        print '\n\n %d set inclusion scores calculated for %s' % (len(J_list_resatm), resatm)
        J_all_resatms[resatm] = J_list_resatm
         
    return J_all_resatms

#############################################################################3#

def extra_set_sim(micros_who_bind_querry, test_resatms = [], test_num = 0):
### micros_who_bind_querry is a dictionary with residue atom types as keys 
### and lists of the line numbers of all that type of residue that are known 
### to bind the fragment as values

    print 'CALCULATING SET INCLUSION FOR NON-BINDING MICRO-ENVIRONMENTS (CONTROL)'

    if test_resatms:
        print 'Checking the requested microenvironment types'
        resatms = [resatm for resatm in test_resatms if resatm in micros_who_bind_querry.keys()]
        print 'Will calculate scores for: \n%s' % (str(resatms))
    else:
        resatms = micros_who_bind_querry.keys()

    # initialize the set inclusion score resatm dictionary
    J_all_resatms = {}

    for resatm in resatms:
        
        # Prepare the data for this resatm type
        print 'PRE-PROCESSING DATA FOR CALCULATION\n'
        frag_binding_ff_dict, non0stdev, resatm_clusters, resatmKB = (
            prep_set_data(resatm, micros_who_bind_querry[resatm]))

        print 'DETERMINING CANDIDATE SET\n'
        # Find the IDs of all the micros not known to bind the fragment
        nonbinding_ID_set = set([micro_ID for micro_ID in range(len(resatmKB)) 
                                  if micro_ID not in set(frag_binding_ff_dict.keys())])
        
        # Determine the candidates (drawn from non-binding set) to compare to the frag binding set
        if test_num != 0:
            # Determine the number of candidates to test
            true_test_num = min(test_num, len(nonbinding_ID_set))
            # Grab a random subset of these micro IDs to be candidates
            candidate_ID_set = random.sample(nonbinding_ID_set, true_test_num)
            # Generate a dictionary of ID# : ff_vector for the candidate microenvironments
            print 'building candidate dictionary'
            candidate_ff_dict = {candidate_ID : resatmKB[candidate_ID] for candidate_ID in candidate_ID_set}
            print 'Will test the inclusion of %d candidates\n' % (len(candidate_ff_dict))
        else:
            candidate_ff_dict = {candidate_ID : resatmKB[candidate_ID] for candidate_ID in nonbinding_ID_set} 
            print 'Will test the inclusion of entire frag-binding set (%d candidates)' % (len(candidate_ff_dict))
        
        # Calculate the inclusion scores of all candidates
        print 'CALCULATING INCLUSION SCORES FOR ALL CANDIDATES\n'
        J_list_resatm = (
            compare_candidates_to_binding_set(candidate_ff_dict, frag_binding_ff_dict, non0stdev, resatm_clusters))
        
        print '\n\n %d set inclusion scores calculated for %s' % (len(J_list_resatm), resatm)
        J_all_resatms[resatm] = J_list_resatm
         
    return J_all_resatms

#############################################################################3#

if __name__ == "__main__":
    import sys
    main(sys.argv[1])
