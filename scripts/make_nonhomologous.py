from argparse import ArgumentParser
from numpy import zeros
import cPickle as pickle
import os

from load_directorypaths import *
from load_datafiles import load_seqclusterid, load_resatm_KB_ann
from load_datafiles import get_cmdline_output
from supportingcode import group_ff_byresatm, make_nonhomologous

###############################################################################
def main():
    parser = ArgumentParser(description="Removes homologous entries from the "
                            "knowledge base using a sequence identity "
                            "threshold defined by --homology")
    parser.add_argument("--folder", type=str, dest="folder",
                        help="folder of input/output files")
    parser.add_argument("--weight", type=float, dest="weight", default=1.0,
                        help="weight to apply to the property standard "
                        "deviation (default: 1.0)")
    parser.add_argument("--homology1", type=int, dest="homology1", default=50,
                        help="max sequence identity permitted between the "
                        "query protein and knowledge base proteins")
    parser.add_argument("--homology2", type=int, dest="homology2", default=50,
                        help="max sequence identity permitted between "
                        "knowledge base proteins")
    parser.add_argument("--pdbfile", type=bool, dest="pdbfile", default=False,
                        help="flag indicating whether the protein is from the "
                        "Protein Data Bank")
    args = parser.parse_args()

    # Read in PDB %ID clusters
    # clustdict[pdbID.chainID] = cluster number
    clustdict1 = load_seqclusterid(str(args.homology1))
    clustdict2 = load_seqclusterid(str(args.homology2))

    # result[i]['ANN']  = for entry i, annotation information
    # result[i]['dVEC'] = for entry i, list of distance to RES.ATM database
    infile = "%s/dissimilaritymatrix.stdev%s.pvar" %(args.folder, args.weight)
    infile = open(infile, 'r')
    result = pickle.load(infile)
    infile.close()

    print ("Making the knowledge base non-homologous at %d%%/%d%% sequence "
           "identity threshold" %(args.homology1, args.homology2))

    # Determine the cluster ID of the protein analyzed (query)
    # Sequence of the query protein
    fastafile = "%s/%s.fasta" %(os.path.dirname(args.folder), 
                                result[0]['ANN'].tag)

    # Blast query protein against the Knowledge Base proteins
    cmd = ("blastp -query %s -db %s/sequences/pdb.seqres.txt -num_alignments 1"
           " -outfmt 6 -evalue 1e-3" %(fastafile, KB_HOME))
    stdout = get_cmdline_output(cmd)[0]
    stdout = stdout.split()
    # Grab the sequence identity (%) of the most similar knowledge base protein
    percentidentity = float(stdout[2])
    if percentidentity > args.homology1:	# A homolog was found
	# Query protein will have the same cluster membership as the homolog
	temp = stdout[1].split("_")
	homolog = "%s.%s" %(temp[0].upper(), temp[1])
	clustid = clustdict1[homolog]
    else:  					# No homolog found
	clustid = -1  # Dummy value

    # Group the feature vectors by RES.ATM tag
    group_byresatm = group_ff_byresatm(result)

    for resatm in group_byresatm:
        # Load annotations for the resatm KB entries 
        anndata = load_resatm_KB_ann(resatm)

        # Determine the clusterID of each database entry
        # Cluster IDs are initialized to -1 for no assignment
        db_clusterids = (zeros((1, len(anndata)), dtype=int) - 1)[0]
        for j in range(0, len(anndata)):
            if anndata[j].tag in clustdict2:
                # If the protein has the same cluster id as the query using 
                # clustdict1, mark this protein for removal (-1)
		if clustdict1[anndata[j].tag] == clustid:
		    db_clusterids[j] = -1
		else:
                    db_clusterids[j] = clustdict2[anndata[j].tag]

        # Find the database entries non-homologous to record i and each other
        for i in group_byresatm[resatm]:
            result[i]['NonH_dVEC'], result[i]['NonH_ANNIDX'] = \
                make_nonhomologous(db_clusterids, result[i]['dVEC'])
            del result[i]['dVEC']

    # Pickle the results
    # result[i]['ANN']         = for entry i, annotation information
    # result[i]['NonH_dVEC']   = for entry i, list of dissimilarity to each 
    # Residue.Atom KB member (non-homologous) (sorted by dissimilarity)
    # result[i]['NonH_ANNIDX'] = for entry i, list of annotation index to 
    # Residue.Atom KB (non-homologous) (sorted by dissimilarity)
    datadir = "%s/ID%d.%d" %(args.folder, args.homology1, args.homology2)
    if not os.path.exists(datadir):
        os.mkdir(datadir)
    outfile = "%s/nonhomologous.pvar" %datadir
    outfile = open(outfile, 'w')
    pickle.dump(result, outfile, protocol=2)
    outfile.close()

    # Remove previous pvar
    infile = "%s/dissimilaritymatrix.stdev%s.pvar" %(args.folder, args.weight)
    #os.remove(infile)


if __name__ == "__main__":
	main()
