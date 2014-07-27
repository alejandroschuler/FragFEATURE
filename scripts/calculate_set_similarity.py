import cPickle as pickle

from load_datafiles import load_resatm_KB_prop, load_resatm_stdev
from supportingcode import dissimilarity_to_KB

def calc_set_sim(micros_who_bind_querry):
#micros_who_bind_querry is a dictionary with residue atom types as keys and lists of the line numbers of all that type of residue that are known to bind the fragment as values

# ASSUMES THAT THE SET DECTIONARY IS LOADED ALREADY

    for resatm in micros_who_bind_querry.keys():
        # Load ff vectors for all micros of the type resatm
        resatmKB = load_resatm_KB_prop(resatm) #I think this imports numpy

        # Remove the vectors for the micros that don't bind our querry fragment
        set_ff = resatmKB[micros_who_bind_querry[resatm],:]

        # Load the pre-calculated standard deviation to use for bit conversion
        # for Residue.Atom
        stdev = load_resatm_stdev(resatm)

	# Keep just the non-zero variance features
        set_ff = set_ff[:,(stdev != 0)]
        non0stdev = stdev[stdev != 0]
        
        for local_micro_num, micro_ff in enumerate(set_ff):
            # Calculate dissimilarity to Residue.Atom knowledge base (KB)
            micro_to_micro_sim_vec = dissimilarity_to_KB(micro_ff, set_ff, non0stdev)
            micro_to_micro_sim_vec = homo_filter(micro_to_micro_sim_vec)
            micro_to_set_sim[resatm][local_micro_num] = set_score(micro_to_micro_sim_vec)

def homo_filter(T)
    

