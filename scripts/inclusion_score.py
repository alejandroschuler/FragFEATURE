import numpy as np
import cPickle as pickle
from load_directorypaths import *
import glob
import pandas as pd
import statsmodels.api as sm
from operator import methodcaller

def main(frag_ID, resatms = []):

    frag_ID = str(frag_ID)
    fragdir = KB_HOME + '/f' + frag_ID   
    if not resatms:
        # Identify the resatm files that are in the drectory
        # so that we don't try to compute on nonexistant data
        file_list = glob.glob(fragdir + '/*tanimotos_full*')
        file_list = map(methodcaller("split", "/"), file_list)
        file_list = [file_name[-1] for file_name in file_list]
        file_info = map(methodcaller("split", "_"), file_list)
        resatms = set([file_datum[0] for file_datum in file_info]) 
    
    # Get the total counts of each resatm to correct the logistic regression
    # for case-control sampling
#    try:
#        resatm_count_dict = pickle.load(open('resatm_count_dict.pvar','w'))
#    except IOError:
#        num_file_name = KB_HOME + '/resatm_counts.txt'
#        num_file = open(num_file_name, 'r')
#        resatm_count_dict = {}
#        for line in num_file:
#            res_data = line.split(" ")
#            count = res_data.pop(0)
#            res_data = res_data[0].split("/")
#            resatm = res_data.pop(0)
#            resatm_count_dict[resatm] = int(count)
#        pickle.dump(resatm_count_dict, open('resatm_count_dict.pvar','w'))
     
    print resatms 
    J_models = dict.fromkeys(resatms)
    J_data = pd.DataFrame(columns =['set_class','score','intercept','resatm'])
    outlist = []
    for resatm in resatms:
        T_intra_file = fragdir + '/' + resatm + '_intra_tanimotos_full.pvar'
        T_extra_file = fragdir + '/' + resatm + '_extra_tanimotos_full.pvar'
        J_models[resatm], J_data_resatm = JCalculator(T_intra_file, T_extra_file)
        J_data_resatm['resatm'] = resatm
        J_data = J_data.append(J_data_resatm, ignore_index=True)

    outfile = fragdir + '/J_data_' + frag_ID + '.csv'
    J_data.to_csv(outfile)
    
    outfile = fragdir + '/J_models_' + frag_ID + '.pvar'
    pickle.dump(J_models, open(outfile,'w'))
                 
def JCalculator(T_intra_file, T_extra_file):
    print '\n\nFrom files:\n%s\n%s\n' % (T_intra_file, T_extra_file)
    T_intra = np.load(T_intra_file)
    T_extra = np.load(T_extra_file)
   
    J_intra = set_score(T_intra)
    J_extra = set_score(T_extra)
    
    n_cases = float(len(J_intra))
    n_controls = float(len(J_extra))

    set_intra = np.ones(n_cases)
    set_extra = np.zeros(n_controls)

    scores = np.concatenate((J_intra, J_extra))
    set_class = np.concatenate((set_intra, set_extra))
    score_table = np.transpose(np.vstack((set_class, scores)))

    J_df = pd.DataFrame(data=score_table, columns = ['set_class', 'score'])
    J_df['intercept'] = 1.0
    J_logit = sm.Logit(J_df['set_class'], J_df[['intercept','score']])          
    J_model = J_logit.fit()

    # Prevalence of the case in the sample
    sample_p = n_cases/(n_cases + n_controls)
    # True prevalence of the case (arbitrary factor set by me. Should change by frag)
    true_p = 1.0/11.0
    # Calculate the corrective factor for the intercept
    correction = np.log(true_p/(1-true_p)) - np.log(sample_p/(1-sample_p))
    # Adjust intercept
    J_model.params.intercept = J_model.params.intercept + correction 
    print 'CORRECTION FACTOR: %f' % (correction) 
    print '\nMODEL WITH CORRECTION:\n'
    print J_model.summary()

    return J_model, J_df

def set_score(T):
    p, k = 100, 0.01
    J = np.sum(T ** p, axis=-1)  ** k
    
    return J

#    def res_summary(self):
#        intra_sm = SummaryObj(self.intra)
#        extra_sm = SummaryObj(self.extra)
#
#        print '\nSummary for Intra-Binding-Set Inclusion Scores'
#        intra_sm.print_summary()
#        print '\nSummary for Extra-Binding-Set Inclusion Scores'
#        extra_sm.print_summary()
#        
#        min_j = min(intra_sm.min, extra_sm.min)
#        max_j = max(intra_sm.max, extra_sm.max)
#        
#        bins = np.linspace(min_j, max_j, 40)
#        intra_hist, bins = np.histogram(self.intra, bins)
#        extra_hist, bins = np.histogram(self.extra, bins)
#        #print '\nbins:\t\t%s' % (str(bins))
#        intra_freq, extra_freq = intra_hist.astype(float)/intra_sm.N, extra_hist.astype(float)/extra_sm.N
#        print '\nHistogram Values'
#        print 'intra:\t%s' % (str(intra_hist))
#        print 'extra:\t%s' % (str(extra_hist))
#
#        B_coef = sum(map(lambda x,y: np.sqrt(x*y), intra_freq, extra_freq))
#        print '\nBhattacharyya Coefficient: %.2f' % (B_coef)
#        B_dist = -np.log(B_coef)
#        print 'Bhattacharyya Distance: %.2f' % (B_dist)
#
#class SummaryObj:
#
#    def __init__(self, vec):
#        self.N = len(vec)
#        self.max = max(vec)
#        self.min = min(vec)
#        self.mean = np.mean(vec)
#        self.std = np.std(vec)
#        self.histogram = np.histogram(vec)
#
#    def print_summary(self):
#        print 'min\t\tmean\t\tmax\t\tstd\t\tsamples'
#        print '%.2f\t\t%.2f\t\t%.2f\t\t%.2f\t\t%d' % (
#            self.min, self.mean, self.max, self.std, self.N)

if __name__ == '__main__':
    import sys
    main(sys.argv[1])
