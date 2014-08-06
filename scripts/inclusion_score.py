import numpy as np
import cPickle as pickle

def main(filename):
    from microenvironment_types import types
    Jobj = JCalculator(filename)
    for resatm in types:
        print '\n-------------------------- %s ------------------------' % (resatm)
        B_dist = Jobj.res_summary(resatm)

class JCalculator:
    
    def __init__(self, T_var):
        if type(T_var) is str:
            # Js will be the name of the file where the vars are picled
            print 'From file: %s' % (T_var)
            T_file = open(T_var, 'r')
            T_dict = pickle.load(T_file)
        else:
            T_dict = T_var   

        self.J = {resatm : dict.fromkeys(['intra','extra']) for resatm in T_dict.keys()}   
        for resatm in T_dict.keys():
             self.J[resatm]['intra'] = self.set_score(T_dict[resatm]['intra'])
             self.J[resatm]['extra'] = self.set_score(T_dict[resatm]['extra'])
                       
    def set_score(self, T):
        p = 5
        k = 4
        J = np.sum(np.abs(T)**p, axis=-1)**k
    
        return J

    def res_summary(self, resatm):
        intra_sm = SummaryObj(self.J[resatm]['intra'])
        extra_sm = SummaryObj(self.J[resatm]['extra'])
       
        print '\nSummary for Intra-Binding-Set Inclusion Scores' 
        intra_sm.print_summary()
        print '\nSummary for Extra-Binding-Set Inclusion Scores' 
        extra_sm.print_summary()
        
        min_j = min(intra_sm.min, extra_sm.min)
        max_j = max(intra_sm.max, extra_sm.max)
        
        N = intra_sm.N
        bins = np.linspace(min_j, max_j, 20)
        intra_hist, bins = np.histogram(self.J[resatm]['intra'], bins)
        extra_hist, bins = np.histogram(self.J[resatm]['extra'], bins)
        #print '\nbins:\t\t%s' % (str(bins))
        print '\nHistogram Values'
        print 'intra:\t%s' % (str(intra_hist))
        print 'extra:\t%s' % (str(extra_hist))

        B_coef = sum(map(lambda x,y: np.sqrt(x*y)/N, intra_hist, extra_hist))
        print '\nBhattacharyya Coefficient: %.2f' % (B_coef)
        B_dist = -np.log(B_coef)
        print 'Bhattacharyya Distance: %.2f' % (B_dist)

        return B_dist

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
    main(sys.argv[0])
