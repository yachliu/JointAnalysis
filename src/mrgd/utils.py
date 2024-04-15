import numba
import numpy as np


@numba.jit(nopython = True)
def pearson_corrcoef_between_arrays(a, b):
    # assert a.shape == b.shape, "Both arrays must have the same shape"
    n, m = a.shape
    corr_matrix = np.zeros((n, n))
    
    for i in range(n):
        for j in range(n):
            mean_ai = np.mean(a[i])
            mean_bj = np.mean(b[j])
            
            cov = np.sum((a[i] - mean_ai) * (b[j] - mean_bj)) / m
            std_ai = np.sqrt(np.sum((a[i] - mean_ai)**2) / m)
            std_bj = np.sqrt(np.sum((b[j] - mean_bj)**2) / m)
            if std_ai == 0 or std_bj == 0:
                corr_matrix[i, j] = 0
            else:
                corr_matrix[i, j] = ((cov / (std_ai * std_bj)) + 1) / 2
                # corr_matrix[i, j] = (cov / (std_ai * std_bj))
    return corr_matrix
