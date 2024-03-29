import numpy as np
import pandas as pd
import scipy as sp
import matplotlib

from statsmodels.nonparametric.kde import KDEUnivariate


def pemp(stat, stat0):
    """ Computes empirical values identically to bioconductor/qvalue empPvals """

    assert len(stat0) > 0
    assert len(stat) > 0

    stat = np.array(stat)
    stat0 = np.array(stat0)

    m = len(stat)
    m0 = len(stat0)

    statc = np.concatenate((stat, stat0))
    v = np.array([True] * m + [False] * m0)
    perm = np.argsort(-statc, kind="mergesort")  # reversed sort, mergesort is stable
    v = v[perm]

    u = np.where(v)[0]
    p = (u - np.arange(m)) / float(m0)

    # ranks can be fractional, we round down to the next integer, ranking returns values starting
    # with 1, not 0:
    ranks = np.floor(sp.stats.rankdata(-stat)).astype(int) - 1
    p = p[ranks]
    p[p <= 1.0 / m0] = 1.0 / m0

    return p

def pi0est(p_values, logger, lambda_ = np.arange(0.05,1.0,0.05), pi0_method = "smoother", smooth_df = 3, smooth_log_pi0 = False):
    """ Estimate pi0 according to bioconductor/qvalue """

    # Compare to bioconductor/qvalue reference implementation
    # import rpy2
    # import rpy2.robjects as robjects
    # from rpy2.robjects import pandas2ri
    # pandas2ri.activate()

    # smoothspline=robjects.r('smooth.spline')
    # predict=robjects.r('predict')

    p = np.array(p_values)

    rm_na = np.isfinite(p)
    p = p[rm_na]
    m = len(p)
    ll = 1
    if isinstance(lambda_, np.ndarray ):
        ll = len(lambda_)
        lambda_ = np.sort(lambda_)

    if (min(p) < 0 or max(p) > 1):
        logger.error("p-values not in valid range [0,1].")
        raise
    elif (ll > 1 and ll < 4):
        logger.errorion("If lambda_ is not predefined (one value), at least four data points are required.")
        raise
    elif (np.min(lambda_) < 0 or np.max(lambda_) >= 1):
        logger.error("Lambda must be within [0,1)")
        raise

    if (ll == 1):
        pi0 = np.mean(p >= lambda_)/(1 - lambda_)
        pi0_lambda = pi0
        pi0 = np.minimum(pi0, 1)
        pi0Smooth = False
    else:
        pi0 = []
        for l in lambda_:
            pi0.append(np.mean(p >= l)/(1 - l))
        pi0_lambda = pi0

        if (pi0_method == "smoother"):
            if smooth_log_pi0:
                pi0 = np.log(pi0)
                spi0 = sp.interpolate.UnivariateSpline(lambda_, pi0, k=smooth_df)
                pi0Smooth = np.exp(spi0(lambda_))
                # spi0 = smoothspline(lambda_, pi0, df = smooth_df) # R reference function
                # pi0Smooth = np.exp(predict(spi0, x = lambda_).rx2('y')) # R reference function
            else:
                spi0 = sp.interpolate.UnivariateSpline(lambda_, pi0, k=smooth_df)
                pi0Smooth = spi0(lambda_)
                # spi0 = smoothspline(lambda_, pi0, df = smooth_df) # R reference function
                # pi0Smooth = predict(spi0, x = lambda_).rx2('y')  # R reference function
            pi0 = np.minimum(pi0Smooth[ll-1],1)
        elif (pi0_method == "bootstrap"):
            minpi0 = np.percentile(pi0,0.1)
            W = []
            for l in lambda_:
                W.append(np.sum(p >= l))
            mse = (np.array(W) / (np.power(m,2) * np.power((1 - lambda_),2))) * (1 - np.array(W) / m) + np.power((pi0 - minpi0),2)
            pi0 = np.minimum(pi0[np.argmin(mse)],1)
            pi0Smooth = False
        else:
            logger.error("pi0_method must be one of 'smoother' or 'bootstrap'.")
            raise
    if (pi0<=0):
        #plot_hist(p, f"p-value density histogram used during pi0 estimation", "p-value", "density histogram", "pi0_estimation_error_pvalue_histogram_plot.pdf")
        logger.error(f"The estimated pi0 <= 0. Check that you have valid p-values or use a different range of lambda. Current lambda range: {lambda_}")
        raise
    return {'pi0': pi0, 'pi0_lambda': pi0_lambda, 'lambda_': lambda_, 'pi0_smooth': pi0Smooth}

def qvalue(p_values, pi0, logger, pfdr = False):
    p = np.array(p_values)

    qvals_out = p
    rm_na = np.isfinite(p)
    p = p[rm_na]

    if (min(p) < 0 or max(p) > 1):
        logger.error("p-values not in valid range [0,1].")
        raise
    elif (pi0 < 0 or pi0 > 1):
        logger.error("pi0 not in valid range [0,1].")
        raise

    m = len(p)
    u = np.argsort(p)
    v = sp.stats.rankdata(p,"max")

    if pfdr:
        qvals = (pi0 * m * p) / (v * (1 - np.power((1 - p), m)))
    else:
        qvals = (pi0 * m * p) / v
    
    qvals[u[m-1]] = np.minimum(qvals[u[m-1]], 1)
    for i in list(reversed(range(0,m-2,1))):
        qvals[u[i]] = np.minimum(qvals[u[i]], qvals[u[i + 1]])

    qvals_out[rm_na] = qvals
    return qvals_out

def bw_nrd0(x, logger):
    if len(x) < 2:
        logger.error("bandwidth estimation requires at least two data points.")
        raise

    hi = np.std(x, ddof=1)
    q75, q25 = np.percentile(x, [75 ,25])
    iqr = q75 - q25
    lo = min(hi, iqr/1.34)
    lo = lo or hi or abs(x[0]) or 1

    return 0.9 * lo *len(x)**-0.2

def lfdr(p_values, pi0, logger, trunc = True, monotone = True, transf = "probit", adj = 1.5, eps = np.power(10.0,-8)):
    """ Estimate local FDR / posterior error probability from p-values according to bioconductor/qvalue """
    p = np.array(p_values)

    # Compare to bioconductor/qvalue reference implementation
    # import rpy2
    # import rpy2.robjects as robjects
    # from rpy2.robjects import pandas2ri
    # pandas2ri.activate()

    # density=robjects.r('density')
    # smoothspline=robjects.r('smooth.spline')
    # predict=robjects.r('predict')

    # Check inputs
    lfdr_out = p
    rm_na = np.isfinite(p)
    p = p[rm_na]

    if (min(p) < 0 or max(p) > 1):
        logger.error("p-values not in valid range [0,1].")
        raise
    elif (pi0 < 0 or pi0 > 1):
        logger.error("pi0 not in valid range [0,1].")
        raise

    # Local FDR method for both probit and logit transformations
    if (transf == "probit"):
        p = np.maximum(p, eps)
        p = np.minimum(p, 1-eps)
        x = sp.stats.norm.ppf(p, loc=0, scale=1)

        # R-like implementation
        bw = bw_nrd0(x, logger)
        myd = KDEUnivariate(x)
        myd.fit(bw=adj*bw, gridsize = 512)
        splinefit = sp.interpolate.splrep(myd.support, myd.density)
        y = sp.interpolate.splev(x, splinefit)
        # myd = density(x, adjust = 1.5) # R reference function
        # mys = smoothspline(x = myd.rx2('x'), y = myd.rx2('y')) # R reference function
        # y = predict(mys, x).rx2('y') # R reference function

        lfdr = pi0 * sp.stats.norm.pdf(x) / y
    elif (transf == "logit"):
        x = np.log((p + eps) / (1 - p + eps))

        # R-like implementation
        bw = bw_nrd0(x, logger)
        myd = KDEUnivariate(x)
        myd.fit(bw=adj*bw, gridsize = 512)

        splinefit = sp.interpolate.splrep(myd.support, myd.density)
        y = sp.interpolate.splev(x, splinefit)
        # myd = density(x, adjust = 1.5) # R reference function
        # mys = smoothspline(x = myd.rx2('x'), y = myd.rx2('y')) # R reference function
        # y = predict(mys, x).rx2('y') # R reference function

        dx = np.exp(x) / np.power((1 + np.exp(x)),2)
        lfdr = (pi0 * dx) / y
    else:
        logger.error("Invalid local FDR method.")
        raise

    if (trunc):
        lfdr[lfdr > 1] = 1
    if (monotone):
        lfdr = lfdr[p.ravel().argsort()]
        for i in range(1,len(x)):
            if (lfdr[i] < lfdr[i - 1]):
                lfdr[i] = lfdr[i - 1]
        lfdr = lfdr[sp.stats.rankdata(p,"min")-1]

    lfdr_out[rm_na] = lfdr
    return lfdr_out

def stats(score_df, score_column, logger):
    score_df_decoy = score_df[score_df["decoy"] == 1]
    stat = score_df[score_column].values
    stat0 = score_df_decoy[score_column].values
    p = pemp(stat, stat0)
    score_df["pvalue"] = p

    pi0 = pi0est(score_df["pvalue"].values, logger)
    qvalues = qvalue(score_df["pvalue"].values, pi0["pi0"], logger)
    pep = lfdr(score_df["pvalue"].values, pi0["pi0"], logger)
    score_df["qvalue"] = qvalues
    score_df["pep"] = pep
    return score_df