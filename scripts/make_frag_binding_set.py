from load_directorypaths import KB_HOME, FragFEATURE_HOME
from microenvironment_types import types
from load_datafiles import load_resatm_KB_boundfrags
import cPickle as pickle
import os

#print KB_HOME
#print FragFEATURE_HOME

def main(querry_frag_string):

    querry_frag_num = int(querry_frag_string)
    micros_who_bind_querry = {}

    for resatm in types:
        print '\nLoading %s binding data...' % (resatm)
        fragdata = load_resatm_KB_boundfrags(resatm)
        print 'Looking for querry-binding microenvironments...'
        microenvironment_list = []
        for microenvironment_num, fraglist  in enumerate(fragdata):
            if querry_frag_num in fraglist:
                microenvironment_list.append(microenvironment_num)
        if microenvironment_list:
            # Only create an entry if there are resatms of the type that bind the fragment
            micros_who_bind_querry[resatm] = microenvironment_list
            print '%d %s fragment-binding microenvironments found' % (
                len(microenvironment_list), resatm)
        else:
            print 'No microenvironments of type %s are known to bind fragment %s' %(resatm, querry_frag_string)

### STORE THE FRAGMENT-BINDING SET ### 

    outdir = '%s/f%s' % (KB_HOME, querry_frag_string)
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    outfile = outdir + '/binding_list.pvar'
    
    print '\n\n Writing results to %s' % (outfile)    
    outfile = open(outfile, 'w')
    pickle.dump(micros_who_bind_querry, outfile)
    outfile.close()

    return

#################################################################################
if __name__ == "__main__":
    import sys
    main(sys.argv[1])

