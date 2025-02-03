from numpy import log, vectorize, delete, array
from collections import defaultdict


def safe_xlogx(s, cutoff=1.0e-08):
    """
    Element-wise function that computes x*log(x) and (correctly) returns zero if x
    is small; 0*log(0) should be interpreted as zero.
    x changed to s
    """
    if s < cutoff:
        return 0.0
    return s * log(s)


def natstobits(s):
    """
    Applies the log change of base formula to convert nat units to bit units.
    """
    return s / log(2)


def entropyd(x, est="ML"):
    """
    Computes the entropy of discrete (integer) data.

    INPUT:
        x: array-like, required
            data for which to compute H[x]

        est: string, optional
            estimator. current options

                'ML' : maximum-likelihood (plugin)
                'MM' : Miller-Maddow corrected
                'JK' : Jackknife estimator (can be slow)

    OUTPUT:
        H[x]: entropy of x, measured in nats
    """
    vec_xlogx = vectorize(safe_xlogx)
    # do the frequency counting
    counts = defaultdict(int)
    for xi in x:
        counts[xi] += 1
    # convert to frequencies
    sumpofx = 1.0 * sum(list(counts.values()))
    pofx = array(list(counts.values())) / sumpofx
    # ML estimator
    H_ML = -1 * vec_xlogx(pofx).sum()
    if est == "ML":
        H = H_ML
    elif est == "MM":
        # nonzero bins have already been removed from pofx
        H = H_ML + (len(pofx) - 1.0) / (2.0 * len(x))
    elif est == "JK":
        Sc = 0
        for i in range(0, len(x)):
            newx = delete(x, i)
            Sc += entropy(newx, bins, "ML")
        H_JK = len(x) * H_ML - ((len(x) - 1.0) / len(x)) * Sc
        H = H_JK
    return H
