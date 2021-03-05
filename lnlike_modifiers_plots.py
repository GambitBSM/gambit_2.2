import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import h5py
from collections import OrderedDict


#
# Read data
#


# Read data for spartan.yaml scan

file_name = "runs/spartan_original/samples/results.hdf5"
group_name = '/spartan/'
f = h5py.File(file_name, 'r')
group = f[group_name]

data_orig = OrderedDict()
data_orig['mu'] = np.array(group['#NormalDist_parameters @NormalDist::primary_parameters::mu'], dtype=np.float64)
data_orig['sigma'] = np.array(group['#NormalDist_parameters @NormalDist::primary_parameters::sigma'], dtype=np.float64)
data_orig['LogLike'] = np.array(group['LogLike'], dtype=np.float64)

LogLike_orig_bf = np.max(data_orig['LogLike'])
print()
print("LogLike_orig_bf :", LogLike_orig_bf)
print("Data set size   :", len(data_orig['LogLike']))


# Read data for spartan.lnlike_modifiers.yaml scan

file_name = "runs/spartan_modified/samples/results.hdf5"
group_name = '/spartan/'

f = h5py.File(file_name, 'r')
group = f[group_name]

data_mod = OrderedDict()
data_mod['mu'] = np.array(group['#NormalDist_parameters @NormalDist::primary_parameters::mu'], dtype=np.float64)
data_mod['sigma'] = np.array(group['#NormalDist_parameters @NormalDist::primary_parameters::sigma'], dtype=np.float64)
data_mod['LogLike'] = np.array(group['LogLike'], dtype=np.float64)
data_mod['ModifiedLogLike'] = np.array(group['ModifiedLogLike'], dtype=np.float64)

LogLike_mod_bf = np.max(data_mod['LogLike'])
print()
print("LogLike_mod_bf :", LogLike_mod_bf)
print("Data set size  :", len(data_mod['LogLike']))
print()


#
# Plot data
#

plt.figure(figsize=(8,12))

# Histograms

nbins=100

plt.subplot(321)

plt.hist2d(data_orig['mu'], data_orig['sigma'], bins=nbins, density=True, norm=mpl.colors.LogNorm())
plt.title('Original\n\nsample density')
plt.xlabel('mu')
plt.ylabel('sigma')


plt.subplot(322)

plt.hist2d(data_mod['mu'], data_mod['sigma'], bins=nbins, density=True, norm=mpl.colors.LogNorm())
plt.title('Modified\n\nsample density')
plt.xlabel('mu')
plt.ylabel('sigma')


# Scatter plots

colormap = plt.get_cmap('viridis')
# bounds = np.arange(-3.1, 0, 0.1)
bounds = [-5.915, -3.09, -1.15, -1e-5, 0]
norm = mpl.colors.BoundaryNorm(bounds, colormap.N)


plt.subplot(323)

ordering = np.argsort(data_orig['LogLike'])
for key in data_orig.keys():
  data_orig[key] = data_orig[key][ordering]

z_data = (data_orig['LogLike'] - data_orig['LogLike'][-1])

plt.scatter(data_orig['mu'], data_orig['sigma'], c=z_data, s=6, cmap=colormap, norm=norm)
plt.title('scanner likelihood')
plt.xlabel('mu')
plt.ylabel('sigma')


plt.subplot(324)

ordering = np.argsort(data_mod['ModifiedLogLike'])
for key in data_mod.keys():
  data_mod[key] = data_mod[key][ordering]

z_data = (data_mod['ModifiedLogLike'] - data_mod['ModifiedLogLike'][-1])

plt.scatter(data_mod['mu'], data_mod['sigma'], c=z_data, s=6, cmap=colormap, norm=norm)
plt.title('scanner likelihood')
plt.xlabel('mu')
plt.ylabel('sigma')



# colormap = plt.get_cmap('viridis')
# bounds = [-9.665, -5.915, -3.09, -1.15, 0]
# bounds = np.arange(-3.1, 0, 0.1)
# norm = mpl.colors.BoundaryNorm(bounds, colormap.N)


plt.subplot(325)

ordering = np.argsort(data_orig['LogLike'])
for key in data_orig.keys():
  data_orig[key] = data_orig[key][ordering]

z_data = (data_orig['LogLike'] - data_orig['LogLike'][-1])

plt.scatter(data_orig['mu'], data_orig['sigma'], c=z_data, s=6, cmap=colormap, norm=norm)
plt.title('physics likelihood')
plt.xlabel('mu')
plt.ylabel('sigma')


plt.subplot(326)

ordering = np.argsort(data_mod['LogLike'])
for key in data_mod.keys():
  data_mod[key] = data_mod[key][ordering]

z_data = (data_mod['LogLike'] - data_mod['LogLike'][-1])

plt.scatter(data_mod['mu'], data_mod['sigma'], c=z_data, s=6, cmap=colormap, norm=norm)
plt.title('physics likelihood')
plt.xlabel('mu')
plt.ylabel('sigma')



plt.tight_layout()
plt.show()
# plt.savefig('lnlike_modifiers.pdf')