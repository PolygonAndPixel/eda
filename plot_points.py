from itertools import product
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.mlab import griddata

import sys

cm = plt.cm.get_cmap('jet')
nDims = 2
three_d = int(sys.argv[3]) > 0
contour = int(sys.argv[4]) > 0
print("Using data from {}".format(sys.argv[1]))
print("Save picture as {}".format(sys.argv[2]))
print("Using 3d projection: {}".format(three_d))
print("Using contour: {}".format(contour))
text = ""
for i in range(5, len(sys.argv)):
    text = text + " " + sys.argv[i]
print("Title of graph is {}".format(text))
base_dir = "output/"
if three_d:
    base_dir += "3D/"
if contour:
    base_dir += "Contour/"

x = []
dim_idx = list(range(2))
dim_names = []
llh = []
with open(sys.argv[1], 'r') as f:
    first = True
    for line in f:
        data = line.split()
        if first:
            dim_names = [str(data[i]) for i in dim_idx]
            first = False                    
        elif data[0] != "Param0":        
            params = [float(data[i]) for i in dim_idx]
            params.append(float(data[nDims]))
            x.append(params)
            llh.append(float(data[nDims]))


x = np.asarray(x)
llh = np.asarray(llh)

if len(x[:,0]) > 3000 and contour:
    step = int((len(x[:,0]) + 3000)/3000)
    x = x[0::step,:]
    llh = llh[0::step]
    
min_llh = np.min(llh)
max_llh = np.max(llh)
levels = np.linspace(min_llh, max_llh, 100)

if three_d:
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    if contour:
        yi = np.linspace(np.min(x[:,0]), np.max(x[:,0]), len(x[:,0])/10)
        zi = np.linspace(np.min(x[:,1]), np.max(x[:,1]), len(x[:,1])/10)
        llh_i = griddata(x[:,0], x[:,1], llh, yi, zi)
        sc = ax.plot_surface(yi, zi, llh_i, linewidth=0, cmap=cm)
    else:
        sc = ax.scatter(x[:,0], x[:,1], llh, 'r', c=llh, cmap=cm, s=3, edgecolors='none')
else:
    fig, ax = plt.subplots()
    if contour:
        yi = np.linspace(np.min(x[:,0]), np.max(x[:,0]), len(x[:,0])/10)
        zi = np.linspace(np.min(x[:,1]), np.max(x[:,1]), len(x[:,1])/10)
        llh_i = griddata(x[:,0], x[:,1], llh, yi, zi)
        sc = plt.contour(yi,zi,llh_i,15,linewidths=0.5,colors='k')
        sc = plt.contourf(yi,zi,llh_i,15,cmap=cm, vmin=min_llh, vmax=max_llh, levels=levels)
        #plt.scatter(x[:,0],x[:,1],marker='o',c='b',s=3)
    else:
        sc = ax.scatter(x[:,0], x[:,1], c=llh, cmap=cm, s=3, edgecolors='none')
    
ax.set_xlabel(dim_names[0])
ax.set_ylabel(dim_names[1])

plt.suptitle(text)
fig.colorbar(sc)
sc.colorbar.set_label('likelihood')
plt.savefig(base_dir + sys.argv[2], dpi=600)
plt.close()

