from argparse import ArgumentParser
from numpy import unique,array,delete
import cPickle as pickle
import itertools

###############################################################################
def main():
    parser = ArgumentParser(description="Define all sets of microenvironments"
                            " that are in spatial proximity (proximity "
                            "indicates within the distance threshold set by "
                            "--dist")
    parser.add_argument("--folder", type=str, dest="folder",
                        help="folder of input/output files")
    parser.add_argument("--dist", type=float, dest="cutoff",
                        help="distance cutoff (angstroms) for defining "
                        "proximal microenvironments")
    parser.add_argument("--setsize", type=int, dest="minsetsize",
                        help="minimum number of microenvironments that must "
                        "be in a set")
    args = parser.parse_args()

    # Read in the microenvironments
    # result[i]['ANN']   = for entry i, annotation information
    # result[i]['FFVEC'] = for entry i, feature vector
    infile = "%s/nonhomologous.pvar" %args.folder
    infile = open(infile, 'r')
    result = pickle.load(infile)
    infile.close()

    # Grab the annotation (FF_Point object) for each microenvironment
    # Grab the similarity of the nearest neighbor for each microenvironment
    entries = []
    nndissimilarity = []
    for i in result:
        entries.append(result[i]['ANN'])
        nndissimilarity.append(round(result[i]['NonH_dVEC'][0]/10000.0, 3))
    N = len(entries)

    # Define all microenvironment sets within the distance cutoff
    neighborslist = []
    for i in range(0, N):
        neighbors = []
        for j in range(0, N):
            dist = entries[i].distance(entries[j])
            if dist <= args.cutoff:
                neighbors.append(j)
        neighborslist.append(neighbors)

    # Generate all subsets of microenvironments
    allsets = []
    for setsize in range(args.minsetsize, N+1):
        combinations = []
        for ffset in neighborslist:
            combinations.extend(list(itertools.combinations(ffset, setsize)))
        combinations = unique(combinations)
        allsets.extend(combinations)

    # Output spatially proximal microenvironment sets and set spread
    outfile = "%s/microenvironment.sets.txt" %args.folder
    outfile = open(outfile, 'w')
    for ffset in allsets:
        # Calculate the max distance between set members
        setsize = len(ffset)
        setspread = -1
	setatoms = ""
        setdissimilarity = []
        for i in range(0,setsize):
            setatoms += ","
            setatoms += entries[ffset[i]].aaInfo
            setdissimilarity.append(nndissimilarity[ffset[i]])
            for j in range(i+1,setsize):
                d = entries[ffset[i]].distance(entries[ffset[j]])
                if d > setspread:
                    setspread = d
        outfile.write("%s\t%.3f\t%s\t%s\n" %(",".join(map(str, ffset)), setspread,
                      setatoms[1:], ",".join(map(str, setdissimilarity))))
    outfile.close()

    print ("Found %d microenvironment sets <- %s angstrom cutoff, minimum set"
           " size of %d" %(len(allsets), args.cutoff, args.minsetsize))


if __name__ == "__main__":
    main()
