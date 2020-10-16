"""
Functions for manipulating evidences
====================================
"""

import numpy as np
from scipy.stats import norm

def combine_runs(*runs):
    """
    @returns Independent runs combined
    """
    mean = np.mean([r[0] for r in runs])
    error = np.sum([r[1]**-2 for r in runs])**-0.5
    return [mean, error]

def add_evidences(*evidences):
    """
    @returns Evidences and errors added together
    """
    sum_ = np.sum([e[0] for e in evidences])
    error = np.sum([e[1]**2 for e in evidences])**0.5
    return [sum_, error]

def bayes_factor(r1, r2):
    """
    @returns Bayes factor
    """
    return np.exp(r1[0] - r2[0])

def z_score(r1, r2):
    """
    @returns Posterior of null converted to Z-score assuming equal prior odds
    """
    return norm.isf(1. / (1. + bayes_factor(r1, r2)))

def evidence(r):
    """
    @returns Evidence from log-evidence
    """
    return np.exp(r[0])
