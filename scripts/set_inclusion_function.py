import numpy  as np
import scipy.misc as sci

def set_score(T):
    b = 1000
    k = 1
    
    result = map(lambda x: b**x, T)
    result = sum(result)
    result = np.log(result)
    result = result - np.log(len(T))
    result = result/np.log(b)
    J = result**k
    # J = (((np.log(sum(map(lambda x: b**x, T)))) - np.log(len(T)))/np.log(b)) ** k 
    
    return J
