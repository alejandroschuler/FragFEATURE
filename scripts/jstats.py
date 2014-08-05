import numpy as np
import cPickle as pickle

def main():
    from microenvironment_types import types
    mfl = 'j_f8027_p5k4.pvar'
    mfl = open(mfl, 'r')
    J = JStats(pickle.load(mfl))
    print 'From file: %s' % (mfl)
    for resatm in types:
        print '\n-------------------------- %s ------------------------' % (resatm)
        B_dist = J.res_summary(resatm)

class JStats:
    
    def __init__(self, Js):
        if type(Js) is str:
            # Js will be the name of the file where the vars are picled
            self.intra, self.extra = pickle.load(Js)
        elif type(Js) is tuple:
            self.intra, self.extra = Js
            
    def res_summary(self, resatm):
        intra_sm = SummaryObj(self.intra[resatm])
        extra_sm = SummaryObj(self.extra[resatm])
       
        print '\nSummary for Intra-Binding-Set Inclusion Scores' 
        intra_sm.print_summary()
        print '\nSummary for Extra-Binding-Set Inclusion Scores' 
        extra_sm.print_summary()
        
        min_j = min(intra_sm.min, extra_sm.min)
        max_j = max(intra_sm.max, extra_sm.max)
        
        N = intra_sm.N
        bins = np.linspace(min_j, max_j, 20)
        intra_hist, bins = np.histogram(self.intra[resatm], bins)
        extra_hist, bins = np.histogram(self.extra[resatm], bins)
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
    main()
