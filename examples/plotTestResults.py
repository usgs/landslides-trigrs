
import h5py
import numpy as np

baseName = "data/outputTRfs_min_ncc"

nSteps = 4
nx = 571
ny = 834

s = (nSteps, nx, ny)

with h5py.File(baseName+".h5", 'r') as f:
    values = np.asarray(f['values'])

values[values == -9999.0] = np.nan


ascValues = np.empty(s)
for i in range(4):
    ascValues[i, :, :] = np.reshape(np.loadtxt(baseName+"_{}.asc".format(i+1), skiprows=6), [nx, ny])

ascValues[ascValues == -9999.0] = np.nan

with h5py.File("test.h5", 'r') as f:
    slic = np.reshape(f['values'], [nx, ny])

slic[slic == -9999.0] = np.nan

import matplotlib as mpl
mpl.use('TkAgg')
import matplotlib.pyplot as plt

plt.figure()
plt.subplot(211)
plt.imshow(ascValues[0, :, :])

plt.subplot(212)
values = np.reshape(values, [nSteps, ny, nx], order='F')
plt.imshow(values[0, :, :])
plt.show()

