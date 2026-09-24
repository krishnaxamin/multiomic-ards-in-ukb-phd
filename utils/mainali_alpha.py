""" Means to implement Mainali et al 2022's alpha co-occurrence metric. """

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr, isinstalled

# import base stuff
base = importr('base')
utils = importr('utils')

# import Mainali packages
if not isinstalled('CooccurrenceAffinity'):
    devtools = importr('devtools')
    devtools.install_github('kpmainali/CooccurrenceAffinity')
mainali = importr('CooccurrenceAffinity')


def mainali_alpha(x, tot_a, tot_b, n):
    # given a 2x2 contingency table reading from left-right, top-bottom as a,b,c,d
    # X = a; mA = a+b; mB = a+c; N = total (a+b+c+d)

    ml_alpha_py = robjects.r['ML.Alpha']
    result = ml_alpha_py(x, robjects.IntVector([tot_a, tot_b, n]), lev=0.95)

    alpha_hat = result.rx2('est')[0]
    pval = result.rx2('pval')[0]

    return {'alpha_hat': alpha_hat, 'pval': pval}
