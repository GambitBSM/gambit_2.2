'''
Get the best fit loglike and chi2
'''

import sys
import os

if len(sys.argv) < 2:
    print "Usage: bestfit_loglikes.py <path_to_samples> [<other_samples>]"
    exit()

totalloglike = 0
totalchisq = 0
for j in range(1,len(sys.argv)):

    path = sys.argv[j]

    # Base file name
    base = '.'.join(os.listdir(path)[0].strip().split('.')[:-1])

    # Find number of output files
    n =  len([name for name in os.listdir(path) if os.path.isfile(os.path.join(path, name)) 
         and len(name.split('.')[-1].split('_')) <= 2
         and "live" not in name
         and "info" not in name
         and "txt_txt" not in name])
    print "Number of files = ", n

    # Worst loglike is -5e10
    bestfit = -5e10

    for i in range(0,n):

        # Filename
        datafile = path + '/' + base + '.txt_' + str(i)
        infofile = path + '/' + base + '.txt_info_' + str(i)

        # If there is only one file sometimes has _0 and sometimes it does not
        if n == 1 and not os.path.exists(datafile):
            datafile = path + '/' + base + '.txt'
            infofile = path + '/' + base + '.txt_info'

        # Get LogLike column for file i
        xe1t_col = 0
        rho_col = 0
        with open(infofile) as finfo:
            for line in finfo:
                if "LogLike" in line:
                    loglike_col = int(line.strip().split()[1][:-1])
                if "#lnL_XENON1T_Anomaly_NuisanceParameters" in line:
                    xe1t_col = int(line.strip().split()[1][:-1])
                if "#lnL_rho0" in line:
                    rho_col = int(line.strip().split()[1][:-1])

        # Read loglike for file and compare to bestfit
        with open(datafile) as fdata:
            for line in fdata:
                loglike = float(line.strip().split()[loglike_col-1])
                if xe1t_col > 0:
                    loglike += float(line.strip().split()[xe1t_col-1])
                if rho_col > 0:
                    loglike += float(line.strip().split()[rho_col-1])
                if loglike > bestfit:
                    bestfit = loglike

    # Chi square
    chisq = -2 * bestfit

    # Total loglike and chisq
    totalloglike += bestfit
    totalchisq += chisq

# Average loglike and chisq
totalloglike = totalloglike / (len(sys.argv)-1)
totalchisq = totalchisq / (len(sys.argv)-1)

# Best fit
print "Best fit log-likelhood for ", base, " is ", totalloglike
print "Best fit chi squared for ", base, " is ", totalchisq
