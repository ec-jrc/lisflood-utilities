import numpy as np


def lmoments(values: np.ndarray) -> np.ndarray:
    """
    Compute first 2 L-moments of dataset (on first axis)
    :param values: N-D array of values
    :return: an estimate of the first two sample L-moments
    """

    nmoments = 3

    # we need to have at least four values in order
    # to make a sample L-moments estimation
    nvalues = values.shape[0]
    
    if nvalues < 4:
        raise ValueError(
            "Insufficient number of values to perform sample L-moments estimation"
        )

    # sort the values into ascending order
    values = np.sort(values, axis=0)

    sums = np.zeros((nmoments, *(values.shape[1:])))

    for i in range(1, nvalues + 1):
        z = i
        term = values[i - 1]
        sums[0] = sums[0] + term
        for j in range(1, nmoments):
            z -= 1
            term = term * z
            sums[j] = sums[j] + term
            
    sums[0] = sums[0] / z

    y = float(nvalues)
    z = float(nvalues)
    
    for j in range(1, nmoments):
        y = y - 1.0
        z = z * y
        sums[j] = sums[j] / z

    k = nmoments
    p0 = -1.0
    for _ in range(2):
        ak = float(k)
        p0 = -p0
        p = p0
        temp = p * sums[0]
        for i in range(1, k):
            ai = i
            p = -p * (ak + ai - 1.0) * (ak - ai) / (ai * ai)
            temp = temp + (p * sums[i])
        sums[k - 1] = temp
        k = k - 1

    lmoments = np.zeros((nmoments, *(values.shape[1:])))
    lmoments[0] = sums[0]
    lmoments[1] = sums[1]
    lmoments[2] = sums[2]

    return lmoments



def lmoments_new(values: np.ndarray) -> np.ndarray:
    """
    Compute first 3 L-moments of dataset (on first axis)
    :param values: N-D array of values
    :return: an estimate of the first two sample L-moments
    """

    nmoments = 3

    # we need to have at least four values in order
    # to make a sample L-moments estimation
    n = values.shape[0]
    
    if n < 4:
        raise ValueError(
            "Insufficient number of values to perform sample L-moments estimation"
        )

    # sort the values into ascending order
    values = np.sort(values, axis=0)

    #Array to store results
    lmoments = np.zeros((nmoments, *(values.shape[1:])))
    
    #1st L-moment
    mean = np.mean(values,axis=0)
    lmoments[0] = mean
    
    #2nd L-moment
    sums = np.zeros((1, *(values.shape[1:])))

    for i in range(1, n + 1):
        
        term = values[i - 1]
        sum_term = (2*i-n-1)*term
        sums = sums + sum_term
        
    lmoments[1] = sums / (n*(n-1))
    
    #3nd L-moment
    sums = np.zeros((1, *(values.shape[1:])))

    for i in range(1, n + 1):
        
        term = values[i - 1]
        sum_term = (n**2 + 3*n - 6*n*i + 6*i**2 - 6*i + 2)*term
        sums = sums + sum_term
        
    lmoments[2] = sums / (n*(n-1)*(n-2))
    

    return lmoments

