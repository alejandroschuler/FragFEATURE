from load_directorypaths import KB_HOME, FragFEATURE_HOME
from microenvironment_types import types
from load_datafiles import load_resatm_KB_boundfrags
import cPickle as pickle

#print KB_HOME
#print FragFEATURE_HOME

querry_frag_num = 8027
micros_who_bind_querry = {}

for resatm in types:
    print resatm
    fragdata = load_resatm_KB_boundfrags(resatm)
    microenvironment_list = []
    for microenvironment_num, fraglist  in enumerate(fragdata):
	if querry_frag_num in fraglist:
            print microenvironment_num
            microenvironment_list.append(microenvironment_num)
    micros_who_bind_querry[resatm] = microenvironment_list

### STORE THE FRAGMENT-BINDING SET ### 

outfile_pvar = 'micros_who_bind_' + str(querry_frag_num) + '.pvar'
outfile_pvar = open(outfile_pvar, 'w')
pickle.dump(micros_who_bind_querry, outfile_pvar)
outfile_pvar.close()

### WRITE TO CSV FILE ###

outfile_txt = 'micros_who_bind_' + str(querry_frag_num) + '.txt'
outfile_txt = open(outfile_txt, 'w')
outfile_txt.write(str(micros_who_bind_querry))
outfile_txt.close()

#note that in these files, a reference to microenv #0 means line #1, etc. Everything starts from 0.

print 'DING!'
