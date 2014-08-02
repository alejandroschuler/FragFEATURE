import cPickle as pickle
import numpy  as np
import scipy.misc as sci
import random
from load_datafiles import load_resatm_KB_prop, load_resatm_stdev
from supportingcode import dissimilarity_to_KB
from load_directorypaths import *

def main(querry_ID_string, test_resatms=[], test_num=0):

    infile = '%s/micros_who_bind_%s.pvar' % (KB_HOME, querry_ID_string)    
    infile = open(infile, 'r')
    micros_who_bind_querry = pickle.load(infile)

    J_intra = intra_set_sim(micros_who_bind_querry, test_resatms, test_num)
    J_extra = extra_set_sim(micros_who_bind_querry, test_resatms, test_num)

    #outfile = '%s/intra_set_inclusion_scores_%s' %(KB_HOME, querry_ID_string)
    #outfile = open(outfile, 'w')
    #pickle.dump(J, outfile)
    
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

def prep_set_data(resatm, frag_binding_ID_set):

        print '\n\n########################### %s ################################\n' % (resatm)

        # Load the clusters for the type
        print 'Loading microenvironment cluster set for %s...' % (resatm)
        resatm_clusters = load_resatm_clusters(resatm)
        print '%d microenvironment cluster mappings loaded.\n' % (len(resatm_clusters))
        
        # Load ff vectors for all micros of the type resatm
        print 'Loading knowledgebase FragFEATURE vectors for %s...' % (resatm)
        resatmKB = load_resatm_KB_prop(resatm)
        print '%d microenvironment FragFEATURE vectors loaded.\n' % (len(resatmKB))
        
        # Load the querry-binding set for the resatm type
        print 'Loading querry fragment binding set for %s...' % (resatm)
        print '%d querry-binding microenvironments loaded.\n' % (len(frag_binding_ID_set))
        
        # Load the pre-calculated standard deviation to use for bit conversion
        stdev = load_resatm_stdev(resatm)
        # Pull the indices of non-zero standard deviation features
        non0stdev_index = [index for index, element in enumerate(stdev) if element != 0]
        # Filter the KB features (columns) for constant features and save the result as a list of lists
        print 'matrix to list happening...\n'
        resatmKB = resatmKB[:,non0stdev_index].tolist()
        # Get the standard deviations of all variant featurs
        non0stdev = stdev[non0stdev_index]
        
        # Get the feature lists for the micros that bind our querry fragment
        print 'Subsetting FragFEATURE knowledgebase for fragment-binding microenvironments...'
        frag_binding_ff_dict = {micro_ID : resatmKB[micro_ID] for micro_ID in frag_binding_ID_set}
        print '%d relevant FragFEATURE vectors retrieved.\n' % (len(frag_binding_ff_dict))

        return (frag_binding_ff_dict, non0stdev, resatm_clusters, resatmKB) 
        

#############################################################################3#


def compare_candidates_to_binding_set(candidate_ff_dict, frag_binding_ff_dict, non0stdev, resatm_clusters):

        # Initialize the set inclusion score vector      
        J_vec = []

        for candidate_ID, candidate_ff in candidate_ff_dict.iteritems():
            print '--------------------- %d ---------------------\n' % (candidate_ID)
            # Calculate dissimilarity to Residue.Atom knowledge base (KB)
            print 'Calculating Tanimoto scores for microenvironment %d...' % (candidate_ID)
            T_dict = tanimoto_dict(candidate_ff, frag_binding_ff_dict, non0stdev)
            # remove the (trivially=1) candidate "self-similarity" score, if the candidate is in the binding set
            try:
                del T_dict[candidate_ID]
            except KeyError:
                pass
            # Filter the vector for homology
            print 'Filtering the scores for homology...'
            T_filtered = homo_filter(T_dict, resatm_clusters)
            print 'Returned %d non-homologous scores.\n' % (len(T_filtered))
            print 'min \t mean \t max'
            print '%.3f \t %.3f \t %.3f\n' % (min(T_filtered), np.mean(T_filtered), max(T_filtered))
            # Calculate J-score (set buddy score)
            print 'Calculating set inclusion score...' 
            J_vec.append(set_score(T_filtered))
            print 'J_%d = %.2f\n' % (candidate_ID, J_vec[-1])

        return J_vec


#############################################################################3#

def set_score(T):
    b = 1000
    k = 1
    
    result = map(lambda x: b**x, T)
    result = sum(result)
    result = np.log(result)
    result = result - np.log(len(T))
    result = result/np.log(b)
    J = result**k
    # J = (((np.log(sum(map(lambda x: b**x, T)))) - np.log(len(T)))/np.log(b)) ** k 
    
    return J

#############################################################################3#

def load_resatm_clusters(resatm):
   
    clusterfile = '%s/%s/%s.micro_clusters.pvar' % (KB_HOME, resatm, resatm)
    clusterfile = open(clusterfile, 'r')
    clusters = pickle.load(clusterfile)

    return clusters

#############################################################################3#

def homo_filter(T_dict, clusters):
# T_dict is a dictionary--  micro_ID : T(m*,micro(ID))
# clusters is a dictionary-- micro_ID : cluster_ID 

    best_score_in_cluster = {}
    singletons = []
    # Break up T_dict into disjoint subsets by homology cluster (and the homeless singletons)
    for micro_ID, score in T_dict.iteritems():
        try:
            cluster_ID = clusters[micro_ID]
            if cluster_ID not in best_score_in_cluster:
                # This is the first score analyzed that is from this cluster, so it is the max so far
                best_score_in_cluster[cluster_ID] = score
#                print '%.3f is the first score in cluster %d' % (score, cluster_ID)
            else:
                # Compare current score to old max, replace old max if new score is better
                best_score_in_cluster[cluster_ID] = max(best_score_in_cluster[cluster_ID], score)
#                print 'A CHALLENGER APPEARS: %.3f is the best score in cluster %d' % (best_score_in_cluster[cluster_ID], cluster_ID)
        except KeyError: 
            singletons.append(score)
#            print 'Micro %d (scored at %.2f) has no cluster' % (micro_ID, score)

    # print '\n\n  #Clusters Represented: %d \t #Singletons: %d \n' % (
    #    len(cluster_subsets), len(singletons))
    # print 'Total Microenvironments: %d' % (len(T_dict))
    T_filtered = best_score_in_cluster.values() + singletons

    return T_filtered

#############################################################################3#

def tanimoto_dict(candidate_ff, frag_binding_ff_dict, non0stdev):
# This is just a wrapper for dis.to.kb, but maps the output: [10000,0] -> [0,1] and returns an ID-keyed dictionary
    
    set_ff_matrix = np.asmatrix([micro_ff for micro_ff in frag_binding_ff_dict.values()])
    # print set_ff_matrix
    T = dissimilarity_to_KB(candidate_ff, set_ff_matrix, non0stdev)
    T = T.astype(float)
    T = (10000 - T)/10000
    print '%d scores calculated.' % (len(T))
    print 'min \t mean \t max'
    print '%.3f \t %.3f \t %.3f' % (min(T), np.mean(T), max(T))
    # make a labeled dictionary that refers to resatm IDs
    T_dict = {frag_binding_ff_dict.keys()[i] : T[i] for i in range(len(T))}

    return T_dict

#############################################################################3#

def homo_filter_test(n, resatm):
    
    import random

    clusters = load_resatm_clusters(resatm)
    test_IDs = random.sample(range(0,max(clusters.keys())),n)
    test_dict = {ID : round(random.random(),3) for ID in test_IDs}

    test_result = homo_filter(test_dict, clusters)

#############################################################################3#

if __name__ == "__main__":
    import sys
    main(sys.argv[1])
