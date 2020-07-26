"""
Log evidences and their errors for solar and DM ALPs from arXiv version 1.
"""

from bayes import combine_runs, add_evidences

# Solar ALP

xe1t_alp = [-21.688189583210804, 2.9673343805847650E-002]
xe1t_alp_r = [-24.378976118255807, 1.8241152697386458E-002]
xe1t_alp_r_wd = [-33.594111494987480, 2.7287332214140699E-002]
xe1t_3h_alp = [-21.414395923422873, 2.7845065229389763E-002]
xe1t_3h_alp_r = [-22.709554198241914, 2.7462782936082987E-002]
xe1t_3h_alp_r_wd = [-32.282489675922214, 3.3312049616034484E-002]
alp_r = [-1.7880299067052574, 1.6627014358437418E-002]
alp_r_wd = [-10.971688775591119, 2.6138961933836878E-002]

# DM ALP

dm_xe1t_alp = [-22.825411868459700, 2.5694833418009769E-002]
dm_alp_r = [-1.7617635459349759, 1.6474666730133890E-002]
dm_alp_r_wd = [-12.490727140791936, 1.8117843570924309E-002]

dm_xe1t_alp_r_1 = [-23.876929302207930, 2.7201848433570080E-002]
dm_xe1t_alp_r_2 = [-23.890639427142432, 2.6941502017563342E-002]
dm_xe1t_alp_r = combine_runs(dm_xe1t_alp_r_1, dm_xe1t_alp_r_2)

dm_xe1t_3h_alp_r_1 = [-22.624137278696601, 2.8187148458826420E-002]
dm_xe1t_3h_alp_r_2 = [-22.608437638087608, 2.8106086087190559E-002]
dm_xe1t_3h_alp_r = combine_runs(dm_xe1t_3h_alp_r_1, dm_xe1t_3h_alp_r_2)

dm_xe1t_alp_r_wd = [-33.944544369856203, 3.5237017736609484E-002]
# same as above but log prior on eta
dm_xe1t_alp_r_wd_logeta = [-33.657722431533024, 3.6599876459608022E-002]

dm_xe1t_3h_alp_1 = [-21.663686916310322, 2.5982589183559394E-002]
dm_xe1t_3h_alp_2 = [-21.639798779173070, 2.5910675497805190E-002]
dm_xe1t_3h_alp_3 = [-21.620496816279534, 2.5898632063533277E-002]
dm_xe1t_3h_alp = combine_runs(dm_xe1t_3h_alp_1, dm_xe1t_3h_alp_2, dm_xe1t_3h_alp_3)

dm_xe1t_3h_alp_r_wd = [-33.160967001605101, 3.1275503779969779E-002]

# Background only

xe1t = [-22.557779416855389, 6.9891095192451185E-003]
xe1t_3h = [-20.946664875526281, 2.1350060516804027E-002]
# Trivial model - no integration hence no error
r = [-4.7979133853e-01, 0]
r_wd = [-1.1300194552e+01, 0]
# Combine evidences trivially
xe1t_r = add_evidences(xe1t, r)
xe1t_3h_r = add_evidences(xe1t_3h, r)
xe1t_r_wd = add_evidences(xe1t, r_wd)
xe1t_3h_r_wd = add_evidences(xe1t_3h, r_wd)
