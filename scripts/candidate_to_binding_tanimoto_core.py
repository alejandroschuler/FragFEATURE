import cPickle as pickle
import numpy  as np
import scipy.misc as sci
from load_datafiles import load_resatm_KB_prop, load_resatm_stdev
from supportingcode import dissimilarity_to_KB
from load_directorypaths import *

#############################################################################3#

def prep_set_data(resatm, frag_binding_IDs):

        print '########################### %s ################################\n' % (resatm)

        # Load the clusters for the type
        print 'Loading microenvironment cluster set for %s...' % (resatm)
        resatm_clusters = load_resatm_clusters(resatm)
        print '%d microenvironment cluster mappings loaded.\n' % (len(resatm_clusters))
        # Initialize the homology filter
        HomoFil = HomologyFilter(resatm_clusters)

        # Load ff vectors for all micros of the type resatm
        print 'Loading knowledgebase FragFEATURE vectors for %s...' % (resatm)
        resatmKB = load_resatm_KB_prop(resatm)
        print '%d microenvironment FragFEATURE vectors loaded.\n' % (len(resatmKB))
        # Load the pre-calculated standard deviation to use for bit conversion
        stdev = load_resatm_stdev(resatm)
        # Pull the indices of non-zero standard deviation features
        non0stdev_index = [index for index, element in enumerate(stdev) if element != 0]
        # Filter the KB features (columns) for constant features and save the result as a list of lists
        resatmKB = resatmKB[:,non0stdev_index]
        # Get the standard deviations of all variant featurs
        non0stdev = stdev[non0stdev_index]
        # Initialize the tanimoto calculator
        TanCalc = TanimotoCalculator(frag_binding_IDs, resatmKB, non0stdev)
        
        return (HomoFil, TanCalc) 
        

#############################################################################3#

def tanimoto_matrix(candidate_IDs, frag_binding_IDs, HomoFil, TanCalc, resatm):

    # Initialize the Tanimoto matrix between the candidates and the binding set members    
    T = [] 
    N = len(candidate_IDs)

    for i, candidate_ID in enumerate(candidate_IDs):
        print '--------------------- %s: %d / %d ---------------------\n' % (resatm, i, N)
        # Calculate dissimilarity to Residue.Atom knowledge base (KB)
        print 'Calculating Tanimoto scores for microenvironment %d...' % (candidate_ID)
        T_dict = TanCalc.score(TanCalc.resatmKB[candidate_ID,:])
        # remove the (trivially=1) candidate "self-similarity" score, if the candidate is in the binding set
        try:
            del T_dict[candidate_ID]
        except KeyError:
            pass
        # Filter the vector for homology
        print 'Filtering the scores for homology...'
        T_filtered = HomoFil.filter(T_dict)
        print 'Returned %d non-homologous scores.\n' % (len(T_filtered))
        print 'min \t mean \t max'
        print '%.3f \t %.3f \t %.3f\n' % (min(T_filtered), np.mean(T_filtered), max(T_filtered))
        
        T.append(T_filtered)

    T = np.asmatrix(T)
    return T

#############################################################################3#

class HomologyFilter:

    def __init__(self, clusters_dict):
        self.clusters_dict = clusters_dict

    def filter(self, T_dict):
        # T_dict is a dictionary--  micro_ID : T(m*,micro(ID))
        # self.clusters is a dictionary-- micro_ID : cluster_ID 
        best_score_in_cluster = {}
        singletons = []
        # Break up T_dict into disjoint subsets by homology cluster (and the homeless singletons)
        for micro_ID, score in T_dict.iteritems():
            try:
                cluster_ID = self.clusters_dict[micro_ID]
                if cluster_ID not in best_score_in_cluster:
                # This is the first score analyzed that is from this cluster, so it is the max so far
                    best_score_in_cluster[cluster_ID] = score
                else:
                # Compare current score to old max, replace old max if new score is better
                    best_score_in_cluster[cluster_ID] = max(best_score_in_cluster[cluster_ID], score)            
            except KeyError: 
                singletons.append(score)
        T_filtered = best_score_in_cluster.values() + singletons

        return T_filtered

#############################################################################3#

class TanimotoCalculator:

    def __init__(self, frag_binding_IDs, resatmKB, non0stdev):
        self.resatmKB = resatmKB
        self.set_ff_matrix = resatmKB[frag_binding_IDs,:]
        self.set_IDs = frag_binding_IDs
        self.non0stdev = non0stdev

    def score(self, candidate_ff):
    # Wrapper for dis.to.kb, but maps the output: [10000,0] -> [0,1] and returns an ID-keyed dictionary
        
        T = dissimilarity_to_KB(candidate_ff, self.set_ff_matrix, self.non0stdev)
        T = T.astype(float)
        T = (10000 - T)/10000
        print '%d scores calculated.' % (len(T))
        print 'min \t mean \t max'
        print '%.3f \t %.3f \t %.3f\n' % (min(T), np.mean(T), max(T))
        # make a labeled dictionary that refers to resatm IDs
        T_dict = {self.set_IDs[i] : T[i] for i in range(len(T))}

        return T_dict

#############################################################################3#

def load_resatm_clusters(resatm):
   
    clusterfile = '%s/%s/%s.micro_clusters.pvar' % (KB_HOME, resatm, resatm)
    clusterfile = open(clusterfile, 'r')
    clusters = pickle.load(clusterfile)

    return clusters

#############################################################################3#

#def homo_filter_test(n, resatm):
#    
#    import random
#
#    clusters = load_resatm_clusters(resatm)
#    test_IDs = random.sample(range(0,max(clusters.keys())),n)
#    test_dict = {ID : round(random.random(),3) for ID in test_IDs}
#
#    test_result = homo_filter(test_dict, clusters)

#############################################################################3#
