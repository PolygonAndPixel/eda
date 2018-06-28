from itertools import product
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
from matplotlib.patches import Rectangle
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.mlab import griddata

import sys

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()

cm = plt.cm.get_cmap('jet_r')
print("Using data from {}".format(sys.argv[1]))
print("Save picture as {}".format(sys.argv[2]))
text = ""
for i in range(3, len(sys.argv)):
    text = text + " " + sys.argv[i]
print("Title of graph is {}".format(text))

base_dir = sys.argv[1]
files = []
for r, d, f in os.walk(base_dir):
    for ff in f:
        files.append(ff)

results = {}
file_names = []
for file_name in files:
	file_names.append(file_name)
    open_this = base_dir + file_name
    with open(open_this, 'r') as f:
		eff = []
		n_llh = []
		max_llh = []
		worst_llh = []
        for line in f:
            data = line.split()
            eff.append(float(data[0]))
			n_llh.append(float(data[1]))
			max_llh.append(float(data[2]))
			worst_llh.append(float(data[3]))
        results[file_name] = {}
		results[file_name]["efficiency"] = eff
		results[file_name]["n_llh_calls"] = n_llh
		results[file_name]["max_llh"] = max_llh
		results[file_name]["worst_llh"] = worst_llh

x = np.asarray(x)
width = 0.4
matplotlib.rcParams.update({'font.size': 20})

for name in results:

	ax2.plot(, speedup1, color='black', marker="x")
# ax2.plot(threads, speedup2, color='black', linestyle='--', marker="o",
# 		 fillstyle='none')
x = np.range(len(results))
rects = ax1.bar(x-width, results[name]["n_llh_calls"], width = width, color=(0.7,0.7,0.8), align='edge', label="1st approach")
#approach1 = mpatches.Patch(color=(0.7,0.7,0.8), label='1st approach')
ax1.set_xlabel('Algorithm')
ax1.set_xticklabels(filenames)
ax1.set_ylabel('Number of function evaluations')
ax2.set_ylabel('Sampling efficiency')

plt.suptitle(text)

plt.savefig(base_dir + sys.argv[2], dpi=600)
plt.close()
