import numpy as np
from sklearn.ensemble import RandomForestClassifier, RandomForestRegressor
from sklearn import metrics
from sklearn import cross_validation

from load_datafiles import load_resatm_KB_prop, load_frag_binding
from load_directorypaths import *
from microenvironment_types import types

def main(frag_ID):
    X,y = prep_data('all', frag_ID)
    model = prep_model(trees=100)
    yh,yt,fpr,tpr,thresh = analyze(X,y,model)

    ROC = np.transpose(np.vstack([fpr,tpr,thresh]))
    np.savetxt('roc_'+frag_ID+'.csv', ROC, delimiter=",")

def prep_data(resatm, frag_ID):

    print 'Prepping data...'
    if resatm == 'all':
        X_dict = {}
        y_dict = {}
        binding_set = load_frag_binding(frag_ID)
        for type in types:
            X_dict[type] = load_resatm_KB_prop(type, filter_flag=0)
            y_dict[type] = np.zeros(len(X_dict[type]))
            y_dict[type][binding_set[type]] = 1
        X = np.concatenate(X_dict.values())
        y = np.concatenate(y_dict.values())
    else:
        X = load_resatm_KB_prop(resatm)
        binding_set = load_frag_binding(frag_ID)[resatm]
        # Make a column with the response variable
        y = np.zeros(len(X))
        y[binding_set] = 1
        #print X.describe()
    return X, y 

def prep_model(trees=10, nfeat='sqrt', type='reg', boot=1):
    
    if type == 'reg':
        clf = RandomForestRegressor(n_estimators=trees, max_features=nfeat, bootstrap=boot, 
        verbose=1, n_jobs = -1)
    elif type == 'class':
        clf = RandomForestClassifier(n_estimators=trees, max_features=nfeat, bootstrap=boot, 
        verbose=1, n_jobs = -1)
    return clf

def analyze(X, y, clf, seed=123):   
    
    np.random.seed(seed)

    # randomly pick a training set
    select = np.random.uniform(0, 1, len(X))
    training_set, test_set = select <= 0.7, select > 0.7

    X_train, X_test = X[training_set,:], X[test_set,:]
    y_train, y_test = y[training_set], y[test_set]
    
    clf.fit(X_train, y_train)

    y_hat = clf.predict(X_test)

    #res = pd.crosstab(y_test, y_hat, rownames=['actual'], colnames=['preds'])
    if str(type(clf)) == "<class 'sklearn.ensemble.forest.RandomForestClassifier'>":
        rec = metrics.recall_score(y_test, y_hat)
        print rec
    elif str(type(clf)) == "<class 'sklearn.ensemble.forest.RandomForestRegressor'>":
        fpr, tpr, thresholds = metrics.roc_curve(y_test, y_hat, pos_label=1)
        print fpr
        print tpr
        print thresholds
    return y_hat, y_test, fpr, tpr, thresholds

if __name__ == '__main__':
    import sys
    main(sys.argv[1])
