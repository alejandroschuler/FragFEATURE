import cPickle as pickle
import numpy  as np
import scipy.misc as sci
from load_datafiles import load_resatm_KB_prop, load_resatm_stdev
from supportingcode import dissimilarity_to_KB
from load_directorypaths import *

def main(querry_ID_string):

    infile = '%s/micros_who_bind_%s.pvar' % (KB_HOME, querry_ID_string)    
    infile = open(infile, 'r')
    micros_who_bind_querry = pickle.load(infile)
    J = within_set_sim(micros_who_bind_querry)

    outfile = '%s/intra_set_inclusion_scores_%s' %(KB_HOME, querry_ID_string)
    outfile = open(outfile, 'w')
    pickle.dump(J, outfile)

#############################################################################3#

def within_set_sim(micros_who_bind_querry):
### micros_who_bind_querry is a dictionary with residue atom types as keys 
### and lists of the line numbers of all that type of residue that are known 
### to bind the fragment as values

    J_all_resatms = {}

    for resatm in micros_who_bind_querry.keys():
        
        # Load the querry-binding set for the resatm type
        print '\n\n################################### %s #####################################\n' % (resatm)
        print 'Loading querry fragment binding set for %s...' % (resatm)
        resatm_ID_set = micros_who_bind_querry[resatm]
        print '%d querry-binding microenvironments loaded.\n' % (len(resatm_ID_set))
        # Load the clusters for the type
        print 'Loading microenvironment cluster set for %s...' % (resatm)
        resatm_clusters = load_resatm_clusters(resatm)
        print '%d microenvironment cluster mappings loaded.\n' % (len(resatm_clusters))
        # Load ff vectors for all micros of the type resatm
        print 'Loading knowledgebase FragFEATURE vectors for %s...' % (resatm)
        resatmKB = load_resatm_KB_prop(resatm)
        print '%d microenvironment FragFEATURE vectors loaded.\n' % (len(resatmKB))
        # Remove the vectors for the micros that don't bind our querry fragment
        # note that now the line number is a local index, to get the global ID, need to do:
        # resatm_ID_set[local_index]
        print 'Subsetting FragFEATURE knowledgebase for fragment-binding microenvironments...'
        set_ff = resatmKB[resatm_ID_set,:]
        print '%d relevant FragFEATURE vectors retrieved.\n' % (len(set_ff))
        # Load the pre-calculated standard deviation to use for bit conversion
        stdev = load_resatm_stdev(resatm)
	    # Keep just the non-zero variance features
        set_ff = set_ff[:,(stdev != 0)]
        non0stdev = stdev[stdev != 0]
        
        J_vec = []

        for candidate_index, candidate_ff in enumerate(set_ff):

            print '---------------------------------- %s : %d : %d -----------------------------------\n' % (resatm, resatm_ID_set[candidate_index], candidate_index)
            # Calculate dissimilarity to Residue.Atom knowledge base (KB)
            print 'Calculating Tanimoto scores for microenvironment %d...' % (candidate_index)
            T_vec = dissimilarity_to_KB(candidate_ff, set_ff, non0stdev)
            print '%d scores calculated.' % (len(T_vec))
            print 'min \t mean \t max'
            print '%.3f \t %.3f \t %.3f' % (min(T_vec), np.mean(T_vec), max(T_vec))
            # make a labeled dictionary that refers to resatm IDs
            T_dict = {resatm_ID_set[i] : T_vec[i] for i in range(len(T_vec))}
            # remove the (trivially=1) candidate "self-similarity" score
            del T_dict[resatm_ID_set[candidate_index]]
            # Filter the vector for homology
            print 'Filtering the scores for homology...'
            T_filtered = homo_filter(T_dict, resatm_clusters)
            print 'Returned %d non-homologous scores.\n' % (len(T_filtered))
            print 'min \t mean \t max'
            print '%.3f \t %.3f \t %.3f' % (min(T_filtered), np.mean(T_filtered), max(T_filtered))
            # Calculate J-score (set buddy score)
            print 'Calculating set inclusion score...' 
            J_vec.append(set_score(T_filtered))
            print 'J_%d = %.2f' % (candidate_index, J_vec[-1])

        print '%d set similarity scores calculated' % (len(J_vec))
        J_all_resatms[resatm] = J_vec

    return J_all_resatms

#############################################################################3#

def set_score(T):
    J = sci.logsumexp(T)
    
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

    cluster_subsets = {}
    singletons = []
    # Break up T_dict into disjoint subsets by homology cluster (and the homeless singletons)
    for micro_ID, score in T_dict.iteritems():
        if micro_ID in clusters.keys():
            cluster_ID = clusters[micro_ID]
            if cluster_ID not in cluster_subsets:
                cluster_subsets[cluster_ID] = [score]
            else:
                cluster_subsets[cluster_ID].append(score)
            # print 'Micro %d (scored at %.2f) added to cluster %d' % (micro_ID, score, cluster_ID)
        else: 
            singletons.append(score)
            # print 'Micro %d (scored at %.2f) has no cluster' % (micro_ID, score)

    # Find the top score in each cluster, toss the rest, since they are not independent votes
    top_score_from_each_cluster = []
    # print '\n\nclust \t best \t scores'
    for cluster_ID, scores in cluster_subsets.iteritems():
        top_score = max(scores)
        top_score_from_each_cluster.append(top_score)
        # print '%d \t %.3f \t %s' % (cluster_ID, top_score, str(scores))
    
    # print '\n\n  #Clusters Represented: %d \t #Singletons: %d \n' % (
    #    len(cluster_subsets), len(singletons))
    # print 'Total Microenvironments: %d' % (len(T_dict))
    T_filtered = top_score_from_each_cluster + singletons

    return T_filtered

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
