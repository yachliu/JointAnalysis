import numba
import numpy as np

@numba.jit(nopython = True)
def rbf_kernel_numba_matrix(x, y, sigma):
    m = x.shape[0]
    n = y.shape[0]
    k = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            distance_squared = (x[i] - y[j]) ** 2
            k[i, j] = np.exp(-distance_squared / (2 * sigma**2))
    return k

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
    for i in range(n):
        corr_matrix[i, i ] = corr_matrix[i, i ] * (n - 1)
    corr_matrix = np.sum(corr_matrix) / ((2 * n - 2) * n)
    return corr_matrix

@numba.jit(nopython = True)
def return_sm_matrix(mr_xics):
    n = mr_xics.shape[0]
    sm_matrix = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            sm_matrix[i, j] = pearson_corrcoef_between_arrays(mr_xics[i], mr_xics[j])
    return sm_matrix

@numba.jit(nopython = True)
def rbf_kernel_numba_matrix(x, y, sigma):
    m = x.shape[0]
    n = y.shape[0]
    k = np.zeros((m, n))
    for i in range(m):
        for j in range(n):
            distance_squared = (x[i] - y[j]) ** 2
            k[i, j] = np.exp(-distance_squared / (2 * sigma**2))
    return k

@numba.jit(nopython = True)
def return_sm_vector(mr_xics, t_xic, t_location):
    n = mr_xics.shape[0]
    sm_vector = np.zeros((n))
    for i in range(n):
        if i != t_location:
            sm_vector[i] = pearson_corrcoef_between_arrays(mr_xics[i], t_xic)
        else:
            sm_vector[i] = pearson_corrcoef_between_arrays(t_xic, t_xic)
    return sm_vector

@numba.jit(nopython = True)
def rbf_kernel_numba_vector(x, y, y_location, sigma):
    n = x.shape[0]
    k = np.zeros((n))
    for i in range(n):
        if i != y_location:
            distance_squared = (x[i] - y) ** 2
            k[i] = np.exp(-distance_squared / (2 * sigma**2))
        else:
            k[i] = 1
    return k

@numba.jit(nopython = True)
def return_reps_vector(repss, y_rep, y_location):
    n = repss.shape[0]
    k = np.empty((n))
    for i in range(n):
        if i != y_location:
            k[i] = repss[i] * y_rep
        else:
            k[i] = y_rep ** 2
    return k