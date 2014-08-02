from load_directorypaths import KB_HOME, FragFEATURE_HOME
from microenvironment_types import types
from load_datafiles import load_resatm_KB_boundfrags
import cPickle as pickle

#print KB_HOME
#print FragFEATURE_HOME

def main(querry_frag_string):

    querry_frag_num = int(querry_frag_string)
    micros_who_bind_querry = {}

    for resatm in types:
        print resatm
        fragdata = load_resatm_KB_boundfrags(resatm)
        microenvironment_list = []
        for microenvironment_num, fraglist  in enumerate(fragdata):
            if querry_frag_num in fraglist:
                print microenvironment_num
                microenvironment_list.append(microenvironment_num)
        if microenvironment_list:
            # Only create an entry if there are resatms of the type that bind the fragment
            micros_who_bind_querry[resatm] = set(microenvironment_list)

### STORE THE FRAGMENT-BINDING SET ### 

    outfile = KB_HOME + '/micros_who_bind_' + str(querry_frag_num) 
    
    outfile_pvar = outfile + '.pvar'
    outfile_pvar = open(outfile_pvar, 'w')
    pickle.dump(micros_who_bind_querry, outfile_pvar)
    outfile_pvar.close()

### WRITE TO CSV FILE ###
#
#    outfile_txt = outfile + '.txt'
#    outfile_txt = open(outfile_txt, 'w')
#    outfile_txt.write(str(micros_who_bind_querry))
#    outfile_txt.close()

#note that in these files, a reference to microenv #0 means line #1, etc. Everything starts from 0.

### this code has been verified manually by looking at lines 0:9 of PHE.PSEU.ALL.boundfrags

    print 'DING!'
    return

#################################################################################
if __name__ == "__main__":
    import sys
    main(sys.argv[1])

