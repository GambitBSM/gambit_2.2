"""
Log evidences and their errors for solar and DM ALPs from arXiv version 2.

We corrected a bug in the inverse Primakoff contributions and changed a prior in the DM ALP
model.

The background-only evidences were not affected. They were repeated anyway though, so I added them.
"""

from bayes import combine_runs, add_evidences
import log_evidences_arxiv_v1

# Solar ALP

xe1t_alp = [-21.572557615201646, 2.9537615986157881E-002]
xe1t_alp_r = [-24.371778884015715, 1.8209124903145329E-002]
xe1t_alp_r_wd = [-33.604092864274854, 2.7188915143322095E-002]
xe1t_3h_alp = [-21.371445105258715, 2.7970153995233428E-002]
xe1t_3h_alp_r = [-22.699052150079545, 2.7360608490005864E-002]
xe1t_3h_alp_r_wd = [-32.197132149906381, 3.3222775591813754E-002]
alp_r = [-1.8237682912342561, 1.6838160775983563E-002]
alp_r = combine_runs(alp_r, log_evidences_arxiv_v1.alp_r)
alp_r_wd = [-10.956788764611288, 2.6094897483242181E-002]
alp_r_wd = combine_runs(alp_r_wd, log_evidences_arxiv_v1.alp_r_wd)

# DM ALP

dm_xe1t_alp = [-22.831325555367659, 2.5555816305341136E-002]
dm_alp_r = [-1.0480841131440342, 1.0955808462168233E-002]
dm_alp_r_wd = [-11.825280645969752, 1.4318926021968518E-002]
dm_xe1t_alp_r = [-23.214541087637048, 2.6374694325378440E-002]
dm_xe1t_3h_alp_r = [-22.077804380214854, 2.5970805132147903E-002]

dm_xe1t_alp_r_wd = [-33.269276958694910, 3.4574068801504043E-002]
# same as above but log prior on eta
dm_xe1t_alp_r_wd_logeta = [-33.035006093991917, 3.5403607941781164E-002]

dm_xe1t_3h_alp = [-21.669943581743322, 2.5987408131907117E-002]
dm_xe1t_3h_alp_r_wd = [-32.633764318944294, 2.9784298143873910E-002]

# Background only

xe1t = [-22.544922904826286, 6.9344717506236112E-003]
xe1t = combine_runs(xe1t, log_evidences_arxiv_v1.xe1t)
xe1t_3h = [-20.896944520740814, 2.1514289552647437E-002]
xe1t_3h = combine_runs(xe1t_3h, log_evidences_arxiv_v1.xe1t_3h)
# Should not have changed (and did not change)
r = log_evidences_arxiv_v1.r
r_wd = log_evidences_arxiv_v1.r_wd

# Combine evidences trivially
xe1t_r = add_evidences(xe1t, r)
xe1t_3h_r = add_evidences(xe1t_3h, r)
xe1t_r_wd = add_evidences(xe1t, r_wd)
xe1t_3h_r_wd = add_evidences(xe1t_3h, r_wd)
