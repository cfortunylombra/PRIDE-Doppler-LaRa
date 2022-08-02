"""
Description: Plot v7v

Author: C. Fortuny-Lombraña
"""
import time
run_time = time.time()
if __name__=="__main__":
    ########################################################################################################################
    ################################################## IMPORT PACKAGES #####################################################
    ########################################################################################################################

    import sys
    sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")
    import os
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib as mpl

    # Import data
    np.set_printoptions(suppress=False,precision=15)
    estimation_information_matrix = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/PODst10_RISEFalse_LaRaTrue_PRIDETrueFalse/estimation_information_matrix.dat"))
    concatenated_link_ends = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/PODst10_RISEFalse_LaRaTrue_PRIDETrueFalse/concatenated_link_ends.dat"))
    concatenated_times = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/PODst10_RISEFalse_LaRaTrue_PRIDETrueFalse/concatenated_times.dat"))

    concatenated_link_ends_long = list()
    concatenated_times_long = list()
    x_partial_long = list()
    y_partial_long = list()
    z_partial_long = list()

    value = 0
    concatenated_link_ends_short = list()
    concatenated_times_short = list()
    x_partial_short = list()
    y_partial_short = list()
    z_partial_short = list()
    for i in range(0,len(concatenated_link_ends)):
        if concatenated_link_ends[i]==value:
            concatenated_link_ends_short.append(concatenated_link_ends[i])
            concatenated_times_short.append(concatenated_times[i])
            x_partial_short.append(estimation_information_matrix[i][8])
            y_partial_short.append(estimation_information_matrix[i][9])
            z_partial_short.append(estimation_information_matrix[i][10])
        else:
            concatenated_link_ends_long.append(concatenated_link_ends_short)
            concatenated_times_long.append(concatenated_times_short)
            x_partial_long.append(x_partial_short)
            y_partial_long.append(y_partial_short)
            z_partial_long.append(z_partial_short)

            concatenated_link_ends_short = list()
            concatenated_times_short = list()
            x_partial_short = list()
            y_partial_short = list()
            z_partial_short = list()
            value = concatenated_link_ends[i]




    plt.figure(figsize=(6,6))
    plt.rcParams.update({'font.size': 12})

    plt.savefig("Partial1.pdf",bbox_inches="tight")
    plt.show()
    plt.close("all")

    '''
    plt.rcParams.update({'font.size': 16})
    # xposition plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    colorline = plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],xposition_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $x$ [m]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xposition_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")
    '''

print("--- %s seconds ---" % (time.time() - run_time))