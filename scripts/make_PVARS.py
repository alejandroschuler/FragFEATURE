from subprocess import Popen
import cPickle as pickle
import os

###############################################################################
folder = os.environ["LOGS"]

# Read in and process resolution values
# resdict[pdbID] = resolution
resfile = folder + "/resolu.idx"
resfile = open(resfile, 'r')
resdict = {}
for line in resfile.readlines()[6:]:  # First 6 lines are header
    line = line.split()
    pdbID = line[0].upper()
    if len(line) == 2:                # No resolution value given
        resdict[pdbID] = 10.0
    else:                             # Resolution value given
        resdict[pdbID] = float(line[2])
resfile.close()

# Pickle the resolution variable
outfile = folder + "/resolu.pvar"
print "Creating %s" %outfile
outfile = open(outfile, 'w')
pickle.dump(resdict, outfile, protocol=2)
outfile.close()

# Process and pickle the protein sequence fasta data
seqfile = folder + "/pdb.seqres.txt"
seqfile = open(seqfile, 'r')
seqdata = seqfile.readlines()
seqfile.close()
fastadict = {}
for i in range(0, len(seqdata), 2):
    pdbID = seqdata[i][1:5].upper()
    chainID = seqdata[i][6]
    tag = "%s.%s" %(pdbID, chainID)
    fastadict[tag] = seqdata[i] + "\n" + seqdata[i+1]

# Pickle the protein sequence variable
outfile = folder + "/pdb.seqres.pvar"
print "Creating %s" %outfile
outfile = open(outfile, 'w')
pickle.dump(fastadict, outfile, protocol=2)
outfile.close()
