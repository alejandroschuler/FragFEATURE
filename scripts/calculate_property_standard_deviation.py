from argparse import ArgumentParser
from numpy import std, shape
import os

from load_directorypaths import *
from load_datafiles import load_resatm_KB_prop

def main():
    ###########################################################################
    parser = ArgumentParser()
    parser.add_argument("--residue", type=str, dest="res",
                        help="target residue")
    parser.add_argument("--atom", type=str, dest="atom",
                        help="target functional atom")
    parser.add_argument("--nr-threshold", type=str, dest="id",
                        help="sequence identity threshold for defining "
                             "'non-redundant'")
    parser.add_argument("--weight", type=float, dest="weight", default=1.0,
                        help="weight to apply to the property standard "
                        "deviation (default: 1.0)")
    args = parser.parse_args()
    ###########################################################################

    # Read in indices corresponding to the non-redundant set
    datadir = "%s/%s.%s/ID%s" %(KB_HOME, args.res, args.atom, args.id)
    infile = "%s/%s.%s.ALL.nr.idx.txt" %(datadir, args.res, args.atom)
    infile = open(infile, 'r')
    nrindices = map(int, infile.read().splitlines())
    infile.close()

    # Load ff vector data (full database) and get the non-redundant set
    ffprop_KB = load_resatm_KB_prop("%s.%s" %(args.res, args.atom))
    ffprop_KB = ffprop_KB[nrindices, :]

    # Calculate the standard deviation of each property
    stdev = std(ffprop_KB, 0)
    for i in range(0, shape(stdev)[1]):
	if stdev[0,i] < 0.0001:  # Handle rounding error
	    stdev[0,i] = 0.0
    stdev = args.weight * stdev

    # Output standard deviation of each property
    stdev = stdev.tolist()[0]
    outfile = "%s/%s.%s.ALL.nr.stdev.%s" %(datadir, args.res, args.atom, 
                                            args.weight)
    outfile = open(outfile,'w')
    outfile.write('\n'.join(map(str,stdev)) + '\n')
    outfile.close()


if __name__ == "__main__":
    main()
