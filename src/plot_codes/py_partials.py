"""
Description: Plot partials

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
    from tudatpy.kernel import constants, numerical_simulation

    # Import data
    np.set_printoptions(suppress=False,precision=15)
    estimation_information_matrix = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/PODst9_RISEFalse_LaRaTrue_PRIDETrueFalse_corr0/estimation_information_matrix.dat"))
    concatenated_link_ends = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/PODst9_RISEFalse_LaRaTrue_PRIDETrueFalse_corr0/concatenated_link_ends.dat"))
    concatenated_times = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/PODst9_RISEFalse_LaRaTrue_PRIDETrueFalse_corr0/concatenated_times.dat"))

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
        if i== len(concatenated_link_ends)-1:
            concatenated_link_ends_short.append(concatenated_link_ends[i])
            concatenated_times_short.append(concatenated_times[i])
            x_partial_short.append(estimation_information_matrix[i][8])
            y_partial_short.append(estimation_information_matrix[i][9])
            z_partial_short.append(estimation_information_matrix[i][10])
            concatenated_link_ends_long.append(concatenated_link_ends_short)
            concatenated_times_long.append(concatenated_times_short)
            x_partial_long.append(x_partial_short)
            y_partial_long.append(y_partial_short)
            z_partial_long.append(z_partial_short)
        elif concatenated_link_ends[i]==value:
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
            concatenated_link_ends_short.append(concatenated_link_ends[i])
            concatenated_times_short.append(concatenated_times[i])
            x_partial_short.append(estimation_information_matrix[i][8])
            y_partial_short.append(estimation_information_matrix[i][9])
            z_partial_short.append(estimation_information_matrix[i][10])

    flat_concatenated_times_long = [item for sublist in concatenated_times_long for item in sublist]

    label_st = dict()
    label_st[0] = "BADARY"
    label_st[1] = "DSS 63"
    label_st[2] = "EFLSBERG"
    label_st[3] = "HARTRAO"
    #label_st[4] = "IRBENE"
    label_st[4] = "MEDICINA"
    label_st[5] = "ONSALA60"
    label_st[6] = "WETTZELL"
    label_st[7] = "WRT0"
    label_st[8] = "YEBES40M"

    plt.figure(figsize=(6,6))
    plt.rcParams.update({'font.size': 16})
    colors = [plt.cm.jet(i) for i in np.linspace(0, 1, len(concatenated_link_ends_long))]
    for i in range(0,len(concatenated_link_ends_long)):
        plt.plot((np.array(concatenated_times_long[i])-min(flat_concatenated_times_long))/constants.JULIAN_DAY,x_partial_long[i],'-o',color=colors[i],linewidth=0.7,markersize=3.5)
    plt.xlabel(r'Time after landing [Earth days]')
    plt.ylabel(r'$\frac{\partial H}{\partial x_{LaRa}}$ [-]')
    plt.grid()
    plt.savefig("Partial1.png",bbox_inches="tight")
    plt.show()
    plt.close("all")

    plt.figure(figsize=(6,6))
    plt.rcParams.update({'font.size': 16})
    colors = [plt.cm.jet(i) for i in np.linspace(0, 1, len(concatenated_link_ends_long))]
    for i in range(0,len(concatenated_link_ends_long)):
        plt.plot((np.array(concatenated_times_long[i])-min(flat_concatenated_times_long))/constants.JULIAN_DAY,y_partial_long[i],'-o',color=colors[i],linewidth=0.7,markersize=3.5)
    plt.xlabel(r'Time after landing [Earth days]')
    plt.ylabel(r'$\frac{\partial H}{\partial y_{LaRa}}$ [-]')
    plt.grid()
    plt.savefig("Partial2.png",bbox_inches="tight")
    plt.show()
    plt.close("all")

    plt.figure(figsize=(6,6))
    plt.rcParams.update({'font.size': 16})
    colors = [plt.cm.jet(i) for i in np.linspace(0, 1, len(concatenated_link_ends_long))]
    for i in range(0,len(concatenated_link_ends_long)):
        plt.plot((np.array(concatenated_times_long[i])-min(flat_concatenated_times_long))/constants.JULIAN_DAY,z_partial_long[i],'-o',label=label_st[int(concatenated_link_ends_long[i][0])],color=colors[i],linewidth=0.7,markersize=3.5)
    plt.legend(loc="upper right")
    plt.xlabel(r'Time after landing [Earth days]')
    plt.ylabel(r'$\frac{\partial H}{\partial z_{LaRa}}$ [-]')
    plt.grid()
    plt.savefig("Partial3.png",bbox_inches="tight")
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