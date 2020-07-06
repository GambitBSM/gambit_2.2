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

# Combine two runs
def combine_two_runs(x,y):
  return np.array([0.5 * (x[0] + y[0]), (x[1]**-2 + y[1]**-2)**-0.5])

dm_xe1t_alp = np.array([-22.825411868459700, 2.5694833418009769E-002])
dm_alp_r_wd = np.array([-12.490727140791936, 1.8117843570924309E-002])

dm_xe1t_alp_r_1 = np.array([-23.876929302207930, 2.7201848433570080E-002])
dm_xe1t_alp_r_2 = np.array([-23.890639427142432, 2.6941502017563342E-002])
dm_xe1t_alp_r = combine_two_runs(dm_xe1t_alp_r_1, dm_xe1t_alp_r_2)

dm_xe1t_3h_alp_r_1 = np.array([-22.624137278696601, 2.8187148458826420E-002])
dm_xe1t_3h_alp_r_2 = np.array([-22.608437638087608, 2.8106086087190559E-002])
dm_xe1t_3h_alp_r  = combine_two_runs(dm_xe1t_3h_alp_r_1, dm_xe1t_3h_alp_r_2) 

dm_xe1t_alp_r_wd = np.array([-33.944544369856203, 3.5237017736609484E-002])

dm_xe1t_3h_alp_1 = np.array([-21.663686916310322,2.5982589183559394E-002])
dm_xe1t_3h_alp_2 = np.array([-21.639798779173070, 2.5910675497805190E-002])
dm_xe1t_3h_alp = combine_two_runs(dm_xe1t_3h_alp_1, dm_xe1t_3h_alp_2)

dm_xe1t_3h_alp_r_wd = np.array([-33.160967001605101, 3.1275503779969779E-002])

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
