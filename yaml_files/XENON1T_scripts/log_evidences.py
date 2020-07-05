"""
Log evidences and their errors for solar and DM ALPs
"""

import numpy as np
from scipy.stats import norm

# Solar ALP

xe1t_alp = np.array([-21.688189583210804, 2.9673343805847650E-002])
xe1t_alp_r = np.array([-24.378976118255807, 1.8241152697386458E-002])
xe1t_alp_r_wd = np.array([-33.594111494987480, 2.7287332214140699E-002])
xe1t_3h_alp = np.array([-21.414395923422873, 2.7845065229389763E-002])
xe1t_3h_alp_r = np.array([-22.709554198241914, 2.7462782936082987E-002])
xe1t_3h_alp_r_wd = np.array([-32.282489675922214, 3.3312049616034484E-002])
alp_r = np.array([-1.7880299067052574, 1.6627014358437418E-002])
alp_r_wd = np.array([-10.971688775591119, 2.6138961933836878E-002])

# DM ALP

dm_xe1t_alp = np.array([-22.825411868459700, 2.5694833418009769E-002])
dm_xe1t_alp_r = np.array([-23.876929302207930, 2.7201848433570080E-002])
dm_xe1t_3h_alp_r = np.array([-22.624137278696601, 2.8187148458826420E-002])
dm_xe1t_alp_r_wd = np.array([-32.706501022154001, 4.0278351498150289E-002])
dm_xe1t_3h_alp = np.array([-21.663686916310322, 2.5982589183559394E-002])

# Combine two runs
dm_xe1t_3h_alp_r_wd_1 = np.array([-32.317710884900158, 3.6813797616004028E-002])
dm_xe1t_3h_alp_r_wd_2 = np.array([-32.305534738325065, 3.6688827942681475E-002])
dm_xe1t_3h_alp_r_wd = np.array([0.5 * (dm_xe1t_3h_alp_r_wd_1[0] + dm_xe1t_3h_alp_r_wd_1[0]),
                                (dm_xe1t_3h_alp_r_wd_1[1]**-2 + dm_xe1t_3h_alp_r_wd_1[1]**-2)**-0.5])

# Background only

xe1t = np.array([-22.557779416855389, 6.9891095192451185E-003])
xe1t_3h = np.array([-20.946664875526281, 2.1350060516804027E-002])
# Trivial model - no integration hence no error
r = np.array([-4.7979133853e-01, 0])
r_wd = np.array([-1.1300194552e+01, 0])
# Combine evidences trivially
xe1t_r = np.array([xe1t[0] + r[0], (xe1t[1]**2 + r[1]**2)**0.5])
xe1t_3h_r = np.array([xe1t_3h[0] + r[0], (xe1t_3h[1]**2 + r[1]**2)**0.5])
xe1t_r_wd = np.array([xe1t[0] + r_wd[0], (xe1t[1]**2 + r_wd[1]**2)**0.5])
xe1t_3h_r_wd = np.array([xe1t_3h[0] + r_wd[0], (xe1t_3h[1]**2 + r_wd[1]**2)**0.5])

# Bayes factors

def bayes_factor(r1, r2):
    return np.exp(r1[0] - r2[0])

def z_score(r1, r2):
    return norm.isf(1. / (1. + bayes_factor(r1, r2)))

def evidence(r):
    return np.exp(r[0])
