import cPickle as pickle
import numpy  as np
from load_datafiles import load_resatm_KB_prop, load_resatm_stdev
from supportingcode import dissimilarity_to_KB
from load_directorypaths import *

def calc_set_sim(micros_who_bind_querry):
#micros_who_bind_querry is a dictionary with residue atom types as keys and lists of the line numbers of all that type of residue that are known to bind the fragment as values

    for resatm in micros_who_bind_querry.keys():
       
        resatm_ID_set = micros_who_bind_querry[resatm] 
        # Load the clusters for the type
        resatm_clusters = load_resatm_clusters(resatm)
        # Load ff vectors for all micros of the type resatm
        resatmKB = load_resatm_KB_prop(resatm)
        # Remove the vectors for the micros that don't bind our querry fragment
        # note that now the line number is a local index, to get the global ID, need to do:
        # resatm_ID_set[local_index]
        set_ff = resatmKB[resatm_ID_set,:]
        # Load the pre-calculated standard deviation to use for bit conversion
        stdev = load_resatm_stdev(resatm)
	    # Keep just the non-zero variance features
        set_ff_matrix = set_ff_matrix[:,(stdev != 0)]
        non0stdev = stdev[stdev != 0]
        
        for candidate_index, micro_ff in enumerate(set_ff):
            # Calculate dissimilarity to Residue.Atom knowledge base (KB)
            T_vec = dissimilarity_to_KB(micro_ff, set_ff, non0stdev)
            # make a labeled dictionary that refers to resatm IDs
            T_dict = {resatm_ID_set[i] : T_vec[i] for i in range(1,len(T_vec))}
            # remove the (trivially=1) candidate "self-similarity" score
            del T_dict[resatm_ID_set[candidate_index]]
            # Filter the vector for homology
            T_filtered = homo_filter(T_dict, resatm_clusters)
            # Calculate J-score (set buddy score)
            J_vec[candidate_index] = set_score(T_filtered)

#############################################################################3#

def set_score(T):
    


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

    cluster_subsets = {}
    singletons = []
    # Break up T_dict into disjoint subsets by homology cluster (and the homeless singletons)
    for micro_ID, score in T_dict.iteritems():
#        try:
#            cluster_ID = clusters[micro_ID]
#            cluster_subsets[cluster_ID].append(score)
#            print 'Micro %d (scored at %.2f) added to cluster %d' % (micro_ID, score, cluster_ID)
#        except KeyError: 
#            singletons.append(score)
#            print 'Micro %d (scored at %.2f) has no cluster' % (micro_ID, score)
        if micro_ID in clusters.keys():
            cluster_ID = clusters[micro_ID]
            if cluster_ID not in cluster_subsets:
                cluster_subsets[cluster_ID] = [score]
            else:
                cluster_subsets[cluster_ID].append(score)
            print 'Micro %d (scored at %.2f) added to cluster %d' % (micro_ID, score, cluster_ID)
        else: 
            singletons.append(score)
            print 'Micro %d (scored at %.2f) has no cluster' % (micro_ID, score)

    # Find the top score in each cluster, toss the rest, since they are not independent votes
    top_score_from_each_cluster = []
    print '\n\nclust \t top \t list'
    for cluster_ID, scores in cluster_subsets.iteritems():
        top_score = max(scores)
        top_score_from_each_cluster.append(top_score)
        print '%d \t %.3f \t %s' % (cluster_ID, top_score, str(scores))
    
    print '\n\n  #Clusters Represented: %d \t #Singletons: %d \n' % (
        len(cluster_subsets), len(singletons))
    print 'Total Microenvironments: %d' % (len(T_dict))
    T_filtered = top_score_from_each_cluster + singletons

    return T_filtered

#############################################################################3#

def homo_filter_test(n, resatm):
    
    import random

    clusters = load_resatm_clusters(resatm)
    test_IDs = random.sample(range(0,max(clusters.keys())),n)
    test_dict = {ID : round(random.random(),3) for ID in test_IDs}

    test_result = homo_filter(test_dict, clusters)

