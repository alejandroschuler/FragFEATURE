import cPickle as pickle
from load_directorypaths import *
from load_datafiles import load_resatm_KB_ann, load_seqclusterid
from microenvironment_types import types
import cPickle as pickle

def main():
    prot_clusters = load_seqclusterid('70') 
    print 'File loaded'
    micro_clusters = cluster_micros(prot_clusters)
    for resatm in types:
        outfile = '%s/%s/%s.micro_clusters.pvar' % (KB_HOME, resatm, resatm)
        outfile = open(outfile,'w')
        pickle.dump(micro_clusters[resatm], outfile)

    return

################################################################################
def invert_dict(dict_by_prot):    

    dict_by_cluster = {}
    for prot, cluster_ID in dict_by_prot.iteritems():
        if cluster_ID not in dict_by_cluster:
            dict_by_cluster[cluster_ID] = [prot]
        else:
            dict_by_cluster[cluster_ID].append(prot)

    return dict_by_cluster

################################################################################
def cluster_micros(dict_by_tag):

    micro_clusters = {}
    for resatm in types:
        micro_to_cluster = {}
        anndata = load_resatm_KB_ann(resatm)
        print "loaded annotations for %s" % (resatm)
        for micro_ID, micro_ann in enumerate(anndata):
            # look up the protein the micro belongs to with micro.tag
            # then look up the cluster that protein belongs to
            try:
                cluster_ID = dict_by_tag[micro_ann.tag] 
                # print "%s %s with tag %s belongs to cluster %s" % (resatm, str(micro_num), micro_ann.tag, str(clusterID))
                micro_to_cluster[micro_ID] = cluster_ID
            except KeyError:
                pass
                # print "%s %s with tag %s has no cluster" % (resatm, str(micro_num), micro_ann.tag)
        # put the dictionary for this micro type into the master dictionary
        micro_clusters[resatm] = micro_to_cluster

    return micro_clusters


################################################################################
if __name__ == "__main__":
    main()
