import numpy as np


def ads(x, a, b, c, d, e, f):
    """asymmetric double sigmoidal (ADS) functionの実装
    Bickerton, G., Paolini, G., Besnard, J. et al. Quantifying the chemical beauty of drugs.
    Nature Chem 4, 90–98 (2012). https://doi.org/10.1038/nchem.1243

    Args:
        x: The value of each CHEM_LABELS
        a-f: asymmetric double sigmoidal (ADS) function parameter
    Returns:
        float: Estimates approximated by ADS
    """
    if d < 0 or e < 0.1 or f < 0.1:
        return np.nan
    exp1 = np.exp(-1 * ((x - c + (d / 2.0)) / e))
    exp2 = np.exp(-1 * ((x - c - (d / 2.0)) / f))
    return a + (((b) / (1.0 + exp1)) * (1.0 - (1.0 / (1.0 + exp2))))
