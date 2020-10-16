"""
A script to output mean deviance and mean parameters from sample chains.
User to specify:
    - /samples directory 
    - column numbers to exclude from mean parameter calculation
    - likelihood column number
"""
import numpy as np 
import argparse
import os
import matplotlib.pyplot as plt
"""
"""
# Import chains from file
def parseArguments():
    # Create argument parser
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--chains_dir", help="At the moment this script takes sample chains as input. Assuming the loglike is the thrid from last column in one of these chains, the DIC is computed.", type=str)
    parser.add_argument("-c", "--exl_columns",nargs='+', type=int, help="Python style list of columns (as defined in the _info_ files) to exclude from mean params calculation. Eg '-c 1 2 3'")    
    parser.add_argument("-L", "--L_column",nargs='+', type=int, help="Integer number of column that likelihood of interest is in.")    
    
    args = parser.parse_args()
    return args
args = parseArguments()
dir_ = args.chains_dir 
chain = []
print("Loading chains...")
for filename in os.listdir(dir_):
    name, file_extension = os.path.splitext(filename)
    if 'info' in file_extension or 'live' in file_extension or '0_txt_0' in file_extension:
        continue
    else:
        print(filename)
        chain.append(np.loadtxt('%s/%s' %(dir_,filename)))
chain =np.vstack(chain)

# Calculate the mean deviance
like_col      = args.L_column[0]
deviance      = -2*(chain[:,like_col-1])
mean_deviance = np.mean(deviance)

# Calculate the mean params (don't include all lnL's for R params etc)
excl  = [x-1 for x in args.exl_columns]
chain = np.delete(chain,np.array(excl),1)
mean_params_o = [np.mean(x) for x in chain.T]

# # Estimates row of chain that gives likelihood at mean params. Not needed if getting likelihoods in gambit.
# params        =  chain.T / np.max( chain.T,axis=0)
# mean_params   = mean_params_o  / np.max( chain.T,axis=0)
# # Estimate the likelihood at the mean params
# def mse(row):
#     res = (row - mean_params)
#     res[np.isnan(res)] = 0
#     return np.dot(res,res)
# msevec = []
# for row in params:
#     msevec.append(mse(row))
# msevec    = np.array(msevec)
# min_val   = min(msevec)
# min_index = np.where(msevec==min_val)
# meanp_est = params_o[min_index]
# print(meanp_est)

# Print mean devance and posterior mean params
print("Mean deviance = ", mean_deviance )
print("Mean params   = ", mean_params_o )