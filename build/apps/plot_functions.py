from itertools import product
import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.mlab import griddata

import sys

def townsend(x):
	t = math.atan2(x[0], x[1])
	constraint = (2*np.cos(t) - 0.5*np.cos(2*t) - 0.25*np.cos(3*t) - 0.125*np.cos(4*t))
	constraint = constraint*constraint
	constraint = constraint + (2*np.sin(t)) * (2*np.sin(t))
	if x[0]*x[0] + x[1]*x[1] < constraint:
		left = np.cos((x[0] - 0.1)*x[1])
		left = left*left
		return -(left - (x[0]*np.sin(3*x[0]+x[1])))
	else:
		return np.inf

def eggholder(x):
	left = -(x[1] + 47)*np.sin(np.sqrt(np.abs(x[0]/2 + (x[1]+47))))
	right = x[0] * np.sin(np.sqrt(np.abs(x[0] - (x[1]+47))))
	return left-right

def rosenbrock(x):
	lh = 0
	lh = lh + (100 * (x[1]-x[0]*x[0]) * (x[1]-x[0]*x[0]) + (x[0]-1)*(x[0]-1))
	return lh

def himmelblau(x):
	left = (x[0]*x[0] + x[1] - 11);
	left = left*left;
	right = (x[0] + x[1]*x[1] - 7);
	right = right*right;
	return left+right

def gauss_shell(x):
	shell_width = 0.1
	r = 2.0
	factor = 1.0/np.sqrt(2*np.pi*shell_width*shell_width);

	left = (x[0]-2.5)*(x[0]-2.5)
	right = (x[0]+2.5)*(x[0]+2.5)
	left = left + x[1]*x[1]
	right = right + x[1]*x[1]

	left = np.sqrt(left) - r
	left = left*left
	left = left/(2 * shell_width*shell_width)
	left = factor * np.exp(-left)

	right = np.sqrt(right) - r
	right = right*right
	right = right/(2 * shell_width*shell_width)
	right = factor * np.exp(-right)

	#if np.isinf(np.log(left + right)):
	#	return 0
	return -(left + right)
	

cm = plt.cm.get_cmap('jet_r')
nDims = 2
three_d = int(sys.argv[2]) > 0
contour = int(sys.argv[3]) > 0
which_func = int(sys.argv[4])
print("Save picture as {}".format(sys.argv[1]))
print("Using 3d projection: {}".format(three_d))
print("Using contour: {}".format(contour))

base_dir = "../../output/"
if three_d:
	base_dir += "3D/"
if contour:
	base_dir += "Contour/"
n_levels = 100
if which_func == 0:
	func = townsend
	down = [-2.25, -2.5]
	up = [2.5001, 1.7501]
	step = 0.01

	file_name = "townsend_"
	title = "Townsend Function"
elif which_func == 1:
	func = eggholder
	down = [-512, -512]
	up = [512.0000, 512.0000]
	step = 1.0
	s = step*0.5
	file_name = "eggholder_"
	title = "Eggholder Function"
elif which_func == 2:
	func = rosenbrock
	down = [-5, -5]
	up = [5.0000, 5.0000]
	step = 0.01
	file_name = "rosenbrock_"
	n_levels = 1000
	title = "Rosenbrock Function"
elif which_func == 3:
	func = himmelblau
	down = [-5, -5]
	up = [5.0000, 5.0000]
	step = 0.01
	file_name = "himmelblau_"
	title = "Himmelblau's Function"
else:
	func = gauss_shell
	down = [-6.0, -6.0]
	up = [6.0000, 6.000]
	step = 0.01
	file_name = "gauss_shell_"
	title = "Gaussian Shells"

x = np.mgrid[down[0]:up[0]:step, down[1]:up[1]:step].reshape(2,-1).T
llh = [func(i) for i in x]

min_llh = np.min(llh)
tmp_llh = np.copy(llh)
mask = np.ma.masked_where(llh == np.finfo(min_llh).max, llh)
max_llh = np.max(mask)
llh = np.asarray(llh)
levels = np.linspace(min_llh, max_llh, n_levels)
if contour:
	x_i = np.arange(down[0],up[0],step)
	y_i = np.arange(down[1],up[1],step)
	llh = llh.reshape(len(x_i), len(y_i) ).T

if three_d:
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	if contour:

		# llh_i = griddata(x[:,0], x[:,1], llh, x[:,0], x[:,1])
		sc = plt.contour(x_i,y_i,llh,15,linewidths=0.5,colors='k')
		sc = plt.contourf(x_i,y_i,llh,15,cmap=cm, vmin=min_llh, vmax=max_llh, levels=levels)
	else:
		sc = ax.scatter(x[:,0], x[:,1], llh, 'r', c=llh, cmap=cm, s=2, edgecolors='none')
	ax.invert_zaxis()
else:
	fig, ax = plt.subplots()
	if contour:
		# yi = np.linspace(np.min(x[:,0]), np.max(x[:,0]), len(x[:,0])/10)
		# zi = np.linspace(np.min(x[:,1]), np.max(x[:,1]), len(x[:,1])/10)
		# llh_i = griddata(x[:,0], x[:,1], llh, x[:,0], x[:,1])
		sc = plt.contour(x_i,y_i,llh,15,linewidths=0.5,colors='k')
		sc = plt.contourf(x_i,y_i,llh,15,cmap=cm, vmin=min_llh, vmax=max_llh, levels=levels)
		# Following line would ploint dots for each point. Use this to show
		# how certain minimizers work by plotting every 100th point or so
		#plt.scatter(x[:,0],x[:,1],marker='o',c='b',s=3)
	else:
		sc = ax.scatter(x[:,0], x[:,1], c=llh, cmap=cm, s=2, edgecolors='none')

ax.set_xlabel("x")
ax.set_ylabel("y")

plt.suptitle(title)
fig.colorbar(sc)
sc.colorbar.set_label('function evaluation')
plt.savefig(base_dir + file_name + sys.argv[1], dpi=600)
plt.close()
