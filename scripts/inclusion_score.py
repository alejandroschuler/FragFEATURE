import numpy as np
import cPickle as pickle
from load_directorypaths import *
import glob

def main(frag_ID):
    from operator import methodcaller

    fragdir = KB_HOME + '/f' + str(frag_ID)    
    file_list = glob.glob(fragdir + '/*tanimotos_full*')
    file_list = map(methodcaller("split", "/"), file_list)
    file_list = [file_name[-1] for file_name in file_list]
    file_info = map(methodcaller("split", "_"), file_list)
    resatms = set([file_datum[0] for file_datum in file_info]) 
   
    print resatms 
    J = dict.fromkeys(resatms)
    outlist = []
    for resatm in resatms:
        T_intra_file = fragdir + '/' + resatm + '_intra_tanimotos_full.pvar'
        T_extra_file = fragdir + '/' + resatm + '_extra_tanimotos_full.pvar'
        J[resatm] = JCalculator(T_intra_file, T_extra_file)
        J[resatm].res_summary()
        for sample in J[resatm].intra:
            outlist.append(resatm + ',intra,' + str(sample))
        for sample in J[resatm].extra:   
            outlist.append(resatm + ',extra,' + str(sample))
        
    outfile = fragdir + '/J_' + frag_ID + '.txt'
    outfile = open(outfile, 'w')
    for item in outlist:
        print>>outfile, item
    
class JCalculator:
    
    def __init__(self, T_intra_file, T_extra_file):
        print 'From files:\n%s\n%s' % (T_intra_file, T_extra_file)
        T_intra = np.load(T_intra_file)
        T_extra = np.load(T_extra_file)
       
        self.intra = self.set_score(T_intra)
        self.extra = self.set_score(T_extra)
                              
    def set_score(self, T):
        p, k = 100, 0.01
        J = np.sum(T ** p, axis=-1)  ** k
        
        return J
                       
    def res_summary(self):
        intra_sm = SummaryObj(self.intra)
        extra_sm = SummaryObj(self.extra)

        print '\nSummary for Intra-Binding-Set Inclusion Scores'
        intra_sm.print_summary()
        print '\nSummary for Extra-Binding-Set Inclusion Scores'
        extra_sm.print_summary()
        
        min_j = min(intra_sm.min, extra_sm.min)
        max_j = max(intra_sm.max, extra_sm.max)
        
        bins = np.linspace(min_j, max_j, 40)
        intra_hist, bins = np.histogram(self.intra, bins)
        extra_hist, bins = np.histogram(self.extra, bins)
        #print '\nbins:\t\t%s' % (str(bins))
        intra_freq, extra_freq = intra_hist.astype(float)/intra_sm.N, extra_hist.astype(float)/extra_sm.N
        print '\nHistogram Values'
        print 'intra:\t%s' % (str(intra_hist))
        print 'extra:\t%s' % (str(extra_hist))

        B_coef = sum(map(lambda x,y: np.sqrt(x*y), intra_freq, extra_freq))
        print '\nBhattacharyya Coefficient: %.2f' % (B_coef)
        B_dist = -np.log(B_coef)
        print 'Bhattacharyya Distance: %.2f' % (B_dist)

class SummaryObj:

    def __init__(self, vec):
        self.N = len(vec)
        self.max = max(vec)
        self.min = min(vec)
        self.mean = np.mean(vec)
        self.std = np.std(vec)
        self.histogram = np.histogram(vec)

    def print_summary(self):
        print 'min\t\tmean\t\tmax\t\tstd\t\tsamples'
        print '%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%d' % (
            self.min, self.mean, self.max, self.std, self.N)

if __name__ == '__main__':
    import sys
    main(sys.argv[1])
