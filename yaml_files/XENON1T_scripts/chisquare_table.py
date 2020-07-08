'''
Compute p-values for solar and dm ALPs
'''

# Background only
xe1t_LL = -22.0015967683
xe1t_r_LL = -22.0015967683 -0.47979133853
xe1t_r_wd_LL = -22.0015967683 -11.300194552

# Background with tritium
xe1t_3h_LL = -17.2404573352
xe1t_3h_r_LL = -17.2404573352-0.47979133853
xe1t_3h_r_wd_LL = -17.2404573352-11.300194552

# Solar ALP
xe1t_alp_LL = -14.7355512286
xe1t_alp_r_LL = -21.7154671554
xe1t_alp_r_wd_LL = -27.75401256

# Solar ALP with tritium
xe1t_3h_alp_LL = -14.7185331839
xe1t_3h_alp_r_LL = -16.8190659867
xe1t_3h_alp_r_wd_LL = -22.8555818742

# DM ALP
dm_xe1t_alp_LL = -13.6026617587
dm_xe1t_alp_r_LL = -13.6310815185
dm_xe1t_alp_r_wd_LL = -21.7347451909

# DM ALP with tritium
dm_xe1t_3h_alp_LL = -12.9617515323
dm_xe1t_3h_alp_r_LL = -13.0176499609
dm_xe1t_3h_alp_r_wd_LL = -21.0347947496

# chisquare
def chisq(x):
    return -2*x

# test statistic defined as tt = -2*(LL(s) - LL(b)), 
# where b is the null hypothesis (bkg or bkg + 3h)
def pval(s,b):
  return -2*(s - b)

tex = ("Background & \\num{{{0}}} & \\num{{{1}}} & \\num{{{2}}} & \\num{{{3}}} & \\num{{{4}}} & \\num{{{5}}} \\\\\n"
    "Background + \ce{{^3H}} & \\num{{{6}}} & \\num{{{7}}} & \\num{{{8}}} & \\num{{{9}}} & \\num{{{10}}} & \\num{{{11}}}\\\\\n"
    "Solar ALP & \\num{{{12}}} & \\num{{{13}}} & \\num{{{14}}} & \\num{{{15}}} & \\num{{{16}}} & \\num{{{17}}}\\\\\n"
    "Solar ALP + \ce{{^3H}} & \\num{{{18}}} & \\num{{{19}}} & \\num{{{20}}} & \\num{{{21}}} & \\num{{{22}}} & \\num{{{23}}}\\\\\n"
    "DM ALP & \\num{{{24}}} & \\num{{{25}}} & \\num{{{26}}} & \\num{{{27}}} & \\num{{{28}}} & \\num{{{29}}}\\\\\n"
    "DM ALP + \ce{{^3H}} & \\num{{{30}}} & \\num{{{31}}} & \\num{{{32}}} & \\num{{{33}}} & \\num{{{34}}} & \\num{{{35}}}\\\\"
    ).format(
    chisq(xe1t_LL), pval(xe1t_LL, xe1t_LL), chisq(xe1t_r_LL), pval(xe1t_r_LL, xe1t_r_LL), chisq(xe1t_r_wd_LL), pval(xe1t_r_wd_LL, xe1t_r_wd_LL),
    chisq(xe1t_3h_LL), pval(xe1t_3h_LL, xe1t_3h_LL), chisq(xe1t_3h_r_LL), pval(xe1t_3h_r_LL, xe1t_3h_r_LL), chisq(xe1t_3h_r_wd_LL), pval(xe1t_3h_r_wd_LL, xe1t_3h_r_wd_LL),
    chisq(xe1t_alp_LL), pval(xe1t_alp_LL, xe1t_LL), chisq(xe1t_alp_r_LL), pval(xe1t_alp_r_LL, xe1t_r_LL), chisq(xe1t_alp_r_wd_LL), pval(xe1t_alp_r_wd_LL, xe1t_r_wd_LL),
    chisq(xe1t_3h_alp_LL), pval(xe1t_3h_alp_LL, xe1t_3h_LL), chisq(xe1t_3h_alp_r_LL), pval(xe1t_3h_alp_r_LL, xe1t_3h_r_LL), chisq(xe1t_3h_alp_r_wd_LL), pval(xe1t_3h_alp_r_wd_LL, xe1t_3h_r_wd_LL),
    chisq(dm_xe1t_alp_LL), pval(dm_xe1t_alp_LL, xe1t_LL), chisq(dm_xe1t_alp_r_LL), pval(dm_xe1t_alp_r_LL, xe1t_r_LL), chisq(dm_xe1t_alp_r_wd_LL), pval(dm_xe1t_alp_r_wd_LL, xe1t_r_wd_LL),
    chisq(dm_xe1t_3h_alp_LL), pval(dm_xe1t_3h_alp_LL, xe1t_3h_LL), chisq(dm_xe1t_3h_alp_r_LL), pval(dm_xe1t_3h_alp_r_LL, xe1t_3h_r_LL), chisq(dm_xe1t_3h_alp_r_wd_LL), pval(dm_xe1t_3h_alp_r_wd_LL, xe1t_3h_r_wd_LL)
    )

print tex
