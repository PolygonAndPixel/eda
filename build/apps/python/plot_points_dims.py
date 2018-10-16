from itertools import product
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import matplotlib.pyplot as plt
import matplotlib
from matplotlib.mlab import griddata

import sys

cm = plt.cm.get_cmap('jet')
nDims = index_t(sys.argv[3])
three_d = False # Old attempt
project = index_t(sys.argv[4]) > 0
contour = index_t(sys.argv[5]) > 0
video = index_t(sys.argv[6]) > 0
min_plot_llh = -1*index_t(sys.argv[7])
print("Using data from {}".format(sys.argv[1]))
print("Save picture as {}".format(sys.argv[2]))
print("Number of dimensions is {}".format(nDims))
print("Using projection: {}".format(project))
print("Using contour: {}".format(contour))
print("Save a video: {}".format(video))
print("Use as minimum llh (only used in projection): {}".format(min_plot_llh))
text = ""
for i in range(8, len(sys.argv)):
    text = text + " " + sys.argv[i]

print("Title of graph is {}".format(text))
if three_d:
    base_dir = "3D/"
else:
    base_dir = ""

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

min_llh = np.min(llh)
max_llh = np.max(llh)
if min_llh > min_plot_llh:
    min_plot_llh = min_llh
levels = np.linspace(min_llh, max_llh, 100)

if video:

    fig, ax = plt.subplots()
    imgs = []
    def update_contour(i, par_pairs, x):
        bins = 4
        parameter_idx = index_t((i/bins)%(len(par_pairs)/2))
        p1, p2 = par_pairs[parameter_idx]

        if p1 == p2 or p1 < p2:
            return 0
        which_bin = i%bins
        other_values = []
        for p3 in dim_idx:
            if p3 != p1 and p3 != p2:
                hist, bin_edges = np.histogram(x[:,p3], bins=bins)
                other_values.append([bin_edges[which_bin], bin_edges[which_bin+1]])
        y = []
        z = []
        llh_2 = []
        other_values = np.asarray(other_values)
        for point in x:
            good = True
            o_idx = 0
            for p3 in dim_idx:
                if p3 != p1 and p3 != p2:
                    if other_values[o_idx, 0] > point[p3] or other_values[o_idx, 1] < point[p3]:
                        good = False
                        break
                    o_idx += 1
            if good:
                llh_2.append(point[nDims])
                y.append(point[p1])
                z.append(point[p2])
        yi = np.linspace(np.min(y), np.max(y), len(y)/10)
        zi = np.linspace(np.min(z), np.max(z), len(z)/10)
        llh_i = griddata((y, z), llh_2, (yi[None,:], zi[:,None]), method='cubic')
        sc = plt.contour(yi,zi,llh_i,15,linewidths=0.5,colors='k')
        sc = plt.contourf(yi,zi,llh_i,15,cmap=plt.cm.jet)
        plt.scatter(y,z,marker='o',c='b',s=5)
        return sc.collections

    # Set up formatting for the movie files
    #Writer = animation.writers['ffmpeg']
    #writer = Writer(fps=15, metadata=dict(artist='Maicon Hieronymus'), bitrate=1800)
    #z = np.array([[1, 1], [1, 1]])
    #sc = plt.contourf([1, 2], [1, 2], z, 15, cmap=plt.cm.jet)
    #anime = animation.FuncAnimation(fig, update_contour, 25,
                                    #fargs=(list(product(dim_idx, dim_idx)), x, sc), interval=50,
                                    #blit=True)
    for i in range(1000):
        te = ax.text(90, 90, text)
        add_this = update_contour(i, list(product(dim_idx, dim_idx)), x)
        if isinstance(add_this, index_t):
            continue
        imgs.append(add_this + [te])
    anime = animation.ArtistAnimation(fig, imgs)
    #plt.show()
    save_dir = base_dir + sys.argv[2]
    FFwriter = animation.FFMpegWriter()
    anime.save('basic_animation.mp4', writer = FFwriter)
    #anime.save(save_dir, writer=writer)

else:
    if project:
        # 6 Parameters make 6*5/2 = 15 plots
        k = 0
        for p1 in dim_idx:

            fig, axs = plt.subplots(ncols=2, nrows=3, sharex='row')
            fig.subplots_adjust(hspace=0.4, wspace=0.4)
            ax_idx = 0
            plt.suptitle(text, fontsize=8)
            plt.tight_layout()
            for p2 in dim_idx:
                print("Plotting for {} and {}".format(dim_names[p1], dim_names[p2]))
                ax = axs[ax_idx//2, ax_idx%2]
                ax_idx += 1

                #y = x[:,p1]
                #z = x[:,p2]
                y = []
                z = []
                llh2 = []
                for i in range(len(llh)):
                    if llh[i] >= min_plot_llh:
                        y.append(x[i,p1])
                        z.append(x[i,p2])
                        llh2.append(llh[i])

                sc = ax.scatter(y, z, marker='o', c=llh2, s=0.5,
                    alpha=0.1, cmap=cm)
                ax.set_xlabel(dim_names[p1], fontsize=6)
                ax.set_ylabel(dim_names[p2], fontsize=6)
                # Make it easier to compare the plots
                ax.set_xlim([min(x[:,p1]), max(x[:,p1])])
                ax.set_ylim([min(x[:,p2]), max(x[:,p2])])
                for tick in ax.xaxis.get_major_ticks():
                    tick.label.set_fontsize(5)
                for tick in ax.yaxis.get_major_ticks():
                    tick.label.set_fontsize(5)

            print(base_dir+sys.argv[2]+ "_" + str(k) +
                "_vmin" + sys.argv[7])
            plt.savefig(base_dir + sys.argv[2] + "_" + str(k) +
                "_vmin" + sys.argv[7], dpi=600)
            k += 1
            plt.close()
        exit()


    for p1, p2 in product(dim_idx, dim_idx):
        if p1 == p2 or p1 < p2:
            continue
        print("Plotting for {} and {}".format(dim_names[p1], dim_names[p2]))

        if project:
            print("Doing nothing")
        #     ax = axs[ax_idx//4, ax_idx%4]
        #     ax_idx += 1
        #
        #     #y = x[:,p1]
        #     #z = x[:,p2]
        #     y = []
        #     z = []
        #     llh2 = []
        #     for i in range(len(llh)):
        #         if llh[i] >= min_plot_llh:
        #             y.append(x[i,p1])
        #             z.append(x[i,p2])
        #             llh2.append(llh[i])
        #
        #     sc = ax.scatter(y, z, marker='o', c=llh2, s=0.5, alpha=0.1, cmap=cm)
        #     ax.set_xlabel(dim_names[p1])
        #     ax.set_ylabel(dim_names[p2])
        #     # Make it easier to compare the plots
        #     ax.set_xlim([min(x[:,p1]), max(x[:,p1])])
        #     ax.set_ylim([min(x[:,p2]), max(x[:,p2])])

            #fig.colorbar(sc)
            #sc.colorbar.set_label('llh')

            #norm = matplotlib.colors.Normalize(vmin=min_llh, vmax=max_llh)
            #cb = matplotlib.colorbar.ColorbarBase(ax, cmap=cm)
            #cb.ax.set_xticklabels(np.arange(min_llh, max_llh, 10))
            # cb = plt.colorbar(sc)
            # cb.set_alpha(1)
            # cb.draw_all()
            # A bug when viewing the plots may be avoided by this
            #cb.solids.set_edgecolor("face")
            #cb.ax.set_xticklabels(np.arange(min_llh, max_llh, 10))

        else:
            bins = 4
            bin_combi = list(range(bins))
            bin_combi = [bin_combi for i in range(len(dim_idx)-2)]
            k = 0
            for all_bins in product(*bin_combi):
                if k < 241:
                    k += 1
                    continue
                if three_d:
                    fig = plt.figure()
                    ax = fig.add_subplot(111, projection='3d')
                else:
                    fig, ax = plt.subplots()
                # Get the most interesting bin for all other dimensions and plot against
                # these fixed values (or values in that range which might not be a good
                # idea if there happens to be a sharp border -> works better with MultiNest
                other_values = []
                b = 0
                for p3 in dim_idx:
                    if p3 != p1 and p3 != p2:
                        idx = all_bins[b]
                        b += 1
                        hist, bin_edges = np.histogram(x[:,p3], bins=bins)
                        other_values.append([bin_edges[idx], bin_edges[idx+1]])
                y = []
                z = []
                llh_2 = []
                other_values = np.asarray(other_values)
                for point in x:
                    good = True
                    o_idx = 0
                    for p3 in dim_idx:
                        if p3 != p1 and p3 != p2:
                            if other_values[o_idx, 0] > point[p3] or other_values[o_idx, 1] < point[p3]:
                                good = False
                                break
                            o_idx += 1
                    if good:
                        llh_2.append(point[nDims])
                        y.append(point[p1])
                        z.append(point[p2])
                if three_d:
                    if contour:
                        yi = np.linspace(np.min(y), np.max(y), len(y)/10)
                        zi = np.linspace(np.min(z), np.max(z), len(z)/10)
                        llh_i = griddata(y, z, llh_2, yi, zi)

                        sc = ax.plot_surface(yi, zi, llh_i, linewidth=0.5, cmap=cm)
                    else:
                        sc = ax.scatter(y, z, llh_2, 'r', c=llh_2, cmap=cm, s=5, edgecolors='none')
                else:
                    if contour:
                        yi = np.linspace(np.min(y), np.max(y), len(y))
                        zi = np.linspace(np.min(z), np.max(z), len(z))
                        llh_i = griddata(y, z, llh_2, yi, zi)
                        sc = plt.contour(yi,zi,llh_i,15,linewidths=0.5,colors='k')
                        sc = plt.contourf(yi,zi,llh_i,15,cmap=cm, vmin=min_llh, vmax=max_llh, levels=levels)
                        plt.scatter(y,z,marker='o',c='b',s=5)

                    else:
                        sc = ax.scatter(y, z, c=llh_2, s=0.5, alpha=0.1, cmap=cm)
                        cb = plt.colorbar(sc)
                        cb.set_alpha(1)
                        cb.draw_all()
            ax.set_xlabel(dim_names[p1])
            ax.set_ylabel(dim_names[p2])
                #for k, vals in enumerate(other_values):
                    #text = text + " " + str(vals[0]) + "-" + str(vals[1])
                    #if k%2 == 0:
                        #text = text + "\n"
            plt.suptitle(text)
            fig.colorbar(sc)
            sc.colorbar.set_label('llh')
            plt.savefig(base_dir + sys.argv[2] + "_" + dim_names[p1] + "_"
                        + dim_names[p2] + "_" + str(k), dpi=600)
            k += 1
            plt.close()
    if three_d:
        plt.show()
