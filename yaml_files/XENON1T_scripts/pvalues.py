'''
Compute p-values for solar and dm ALPs
'''

# Background only
xe1t_LL = -21.962808175
xe1t_r_LL = -21.962808175 -0.47979133853
xe1t_r_wd_LL = -21.962808175 -11.300194552

# Background with tritium
xe1t_3h_LL = -16.582389903
xe1t_3h_r_LL = -16.582389903-0.47979133853
xe1t_3h_r_wd_LL = -16.582389903-11.300194552

# Solar ALP
xe1t_alp_LL = -13.814222014
xe1t_alp_r_LL = -21.706297028
xe1t_alp_r_wd_LL = -27.746468275

# Solar ALP with tritium
xe1t_3h_alp_LL = -13.822769243
xe1t_3h_alp_r_LL = -15.714566705
xe1t_3h_alp_r_wd_LL = -21.74603656

# DM ALP
dm_xe1t_alp_LL = -14.67700751
dm_xe1t_alp_r_LL = -14.6774453905
dm_xe1t_alp_r_wd_LL = -22.784802207

# DM ALP with tritium
dm_xe1t_3h_alp_LL = -13.657778231
dm_xe1t_3h_alp_r_LL = -13.6668861425
dm_xe1t_3h_alp_r_wd_LL = -21.700447156

# chisquare
def chisq(x):
    return -2*x

# test statistic defined as tt = -2*(LL(s) - LL(b)), 
# where b is the null hypothesis (bkg + 3h)
def pval(s,b):
  return -2*(s - b)

tex = ("Background & \\num{{{0}}} & \\num{{{1}}} & \\num{{{2}}} & \\num{{{3}}} & \\num{{{4}}} & \\num{{{5}}} \\\\\n"
    "Background + \ce{{^3H}} & \\num{{{6}}} & \\num{{{7}}} & \\num{{{8}}} & \\num{{{9}}} & \\num{{{10}}} & \\num{{{11}}}\\\\\n"
    "Solar ALP & \\num{{{12}}} & \\num{{{13}}} & \\num{{{14}}} & \\num{{{15}}} & \\num{{{16}}} & \\num{{{17}}}\\\\\n"
    "Solar ALP + \ce{{^3H}} & \\num{{{18}}} & \\num{{{19}}} & \\num{{{20}}} & \\num{{{21}}} & \\num{{{22}}} & \\num{{{23}}}\\\\\n"
    "DM ALP & \\num{{{24}}} & \\num{{{25}}} & \\num{{{26}}} & \\num{{{27}}} & \\num{{{28}}} & \\num{{{29}}}\\\\\n"
    "DM ALP + \ce{{^3H}} & \\num{{{30}}} & \\num{{{31}}} & \\num{{{32}}} & \\num{{{33}}} & \\num{{{34}}} & \\num{{{35}}}\\\\"
    ).format(
    chisq(xe1t_LL), pval(xe1t_LL, xe1t_3h_LL), chisq(xe1t_r_LL), pval(xe1t_r_LL, xe1t_3h_r_LL), chisq(xe1t_r_wd_LL), pval(xe1t_r_wd_LL, xe1t_3h_r_wd_LL),
    chisq(xe1t_3h_LL), pval(xe1t_3h_LL, xe1t_3h_LL), chisq(xe1t_3h_r_LL), pval(xe1t_3h_r_LL, xe1t_3h_r_LL), chisq(xe1t_3h_r_wd_LL), pval(xe1t_3h_r_wd_LL, xe1t_3h_r_wd_LL),
    chisq(xe1t_alp_LL), pval(xe1t_alp_LL, xe1t_3h_LL), chisq(xe1t_alp_r_LL), pval(xe1t_alp_r_LL, xe1t_3h_r_LL), chisq(xe1t_alp_r_wd_LL), pval(xe1t_alp_r_wd_LL, xe1t_3h_r_wd_LL),
    chisq(xe1t_3h_alp_LL), pval(xe1t_3h_alp_LL, xe1t_3h_LL), chisq(xe1t_3h_alp_r_LL), pval(xe1t_3h_alp_r_LL, xe1t_3h_r_LL), chisq(xe1t_3h_alp_r_wd_LL), pval(xe1t_3h_alp_r_wd_LL, xe1t_3h_r_wd_LL),
    chisq(dm_xe1t_alp_LL), pval(dm_xe1t_alp_LL, xe1t_3h_LL), chisq(dm_xe1t_alp_r_LL), pval(dm_xe1t_alp_r_LL, xe1t_3h_r_LL), chisq(dm_xe1t_alp_r_wd_LL), pval(dm_xe1t_alp_r_wd_LL, xe1t_3h_r_wd_LL),
    chisq(dm_xe1t_3h_alp_LL), pval(dm_xe1t_3h_alp_LL, xe1t_3h_LL), chisq(dm_xe1t_3h_alp_r_LL), pval(dm_xe1t_3h_alp_r_LL, xe1t_3h_r_LL), chisq(dm_xe1t_3h_alp_r_wd_LL), pval(dm_xe1t_3h_alp_r_wd_LL, xe1t_3h_r_wd_LL)
    )

print tex
