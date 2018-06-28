# Plot an llh scan like the one by Martin Leuermann
# Set all parameters except for the time to the true value and
# plot the llh with different timings
import numpy as np
import matplotlib.pyplot as plt

import sys

cm = plt.cm.get_cmap('jet')
nDims = int(sys.argv[3])

print("Using data from {}".format(sys.argv[1]))
print("Save picture as {}".format(sys.argv[2]))
print("Number of dimensions is {}".format(nDims))
text = ""
for i in range(8, len(sys.argv)):
    text = text + " " + sys.argv[i]

print("Title of graph is {}".format(text))

x = []
dim_idx = list(range(nDims))
dim_names = []
llh = []
with open(sys.argv[1], 'r') as f:
    for line in f:
        data = line.split()
        try:
            params = [float(data[i]) for i in dim_idx]
            params.append(float(data[nDims]))
            x.append(params)
            llh.append(float(data[nDims]))
        except:
            dim_names = [str(data[i]) for i in dim_idx]


x = np.asarray(x)
max_llh = np.max(llh)
x_true = x[np.argmin(llh)]
print("Names of the parameter:")
print(dim_names)
x_true_list = []
llh_true_list = []
time_idx = 0
for i in range(len(dim_names)):
    if dim_names[i] == "T":
        time_idx = i
        break

# for p, l in zip(x, llh):
#     eps = 1e1
#     close = True
#     for i in range(len(p)):
#         if i != time_idx and not np.allclose([p[i]], [x_true[i]], eps):
#             close = False
#     if close:
#         x_true_list.append(p[time_idx])
#         llh_true_list.append(l)
print(llh_true_list)
plt.figure(2, (16,12))

plt.scatter(x[:,time_idx], llh, label="llh scan")
plt.vlines(0, np.min(llh), max_llh,
    color="b", linewidth=3.0, linestyle="--", label="truth")

# plt.scatter(x_true_list, llh_true_list, label="llh scan")
# plt.vlines(0, np.min(llh_true_list), np.max(llh_true_list),
#     color="b", linewidth=3.0, linestyle="--", label="truth")
plt.legend()
plt.xlim([ -100.0,100.0])
plt.show()
