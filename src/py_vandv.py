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
    xposition = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xposition_vandv_new1.dat"))
    yposition = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/yposition_vandv_new1.dat"))
    zposition = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/zposition_vandv_new1.dat"))
    xdotvelocity = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xdotvelocity_vandv_new1.dat"))
    ydotvelocity = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ydotvelocity_vandv_new1.dat"))
    zdotvelocity = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/zdotvelocity_vandv_new1.dat"))
    corefactor = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/corefactor_vandv_new1.dat"))
    sigmaFCN = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/sigmaFCN_vandv_new1.dat"))
    xLaRa = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xLaRa_vandv_new1.dat"))
    yLaRa = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/yLaRa_vandv_new1.dat"))
    zLaRa = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/zLaRa_vandv_new1.dat"))
    cos1spin = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/cos1spin_vandv_new1.dat"))
    sin1spin = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/sin1spin_vandv_new1.dat"))
    cos2spin = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/cos2spin_vandv_new1.dat"))
    sin2spin = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/sin2spin_vandv_new1.dat"))
    cos3spin = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/cos3spin_vandv_new1.dat"))
    sin3spin = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/sin3spin_vandv_new1.dat"))
    cos4spin = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/cos4spin_vandv_new1.dat"))
    sin4spin = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/sin4spin_vandv_new1.dat"))
    xpcos1 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpcos1_vandv_new1.dat"))
    xpsin1 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpsin1_vandv_new1.dat"))
    ypcos1 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypcos1_vandv_new1.dat"))
    ypsin1 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypsin1_vandv_new1.dat"))
    xpcos2 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpcos2_vandv_new1.dat"))
    xpsin2 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpsin2_vandv_new1.dat"))
    ypcos2 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypcos2_vandv_new1.dat"))
    ypsin2 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypsin2_vandv_new1.dat"))
    xpcos3 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpcos3_vandv_new1.dat"))
    xpsin3 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpsin3_vandv_new1.dat"))
    ypcos3 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypcos3_vandv_new1.dat"))
    ypsin3 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypsin3_vandv_new1.dat"))
    xpcos4 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpcos4_vandv_new1.dat"))
    xpsin4 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpsin4_vandv_new1.dat"))
    ypcos4 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypcos4_vandv_new1.dat"))
    ypsin4 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypsin4_vandv_new1.dat"))
    xpcos5 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpcos5_vandv_new1.dat"))
    xpsin5 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpsin5_vandv_new1.dat"))
    ypcos5 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypcos5_vandv_new1.dat"))
    ypsin5 = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypsin5_vandv_new1.dat"))
    
    # Initialize empty lists
    corr_list = list()
    station_number_list = list()

    xposition_list = list()
    yposition_list = list()
    zposition_list = list()
    xdotvelocity_list = list()
    ydotvelocity_list = list()
    zdotvelocity_list = list()
    corefactor_list = list()
    sigmaFCN_list = list()
    xLaRa_list = list()
    yLaRa_list = list()
    zLaRa_list = list()
    cos1spin_list = list()
    sin1spin_list = list()
    cos2spin_list = list()
    sin2spin_list = list()
    cos3spin_list = list()
    sin3spin_list = list()
    cos4spin_list = list()
    sin4spin_list = list()
    xpcos1_list = list()
    xpsin1_list = list()
    ypcos1_list = list()
    ypsin1_list = list()
    xpcos2_list = list()
    xpsin2_list = list()
    ypcos2_list = list()
    ypsin2_list = list()
    xpcos3_list = list()
    xpsin3_list = list()
    ypcos3_list = list()
    ypsin3_list = list()
    xpcos4_list = list()
    xpsin4_list = list()
    ypcos4_list = list()
    ypsin4_list = list()
    xpcos5_list = list()
    xpsin5_list = list() 
    ypcos5_list = list()
    ypsin5_list = list()

    i = 0
    for station in np.linspace(1,10,10):
        corr_list.append(np.linspace(0,9,10)/10)
        station_number_list.append(station*np.ones(10))

        xposition_list.append(xposition[i:i+10])
        yposition_list.append(yposition[i:i+10])
        zposition_list.append(zposition[i:i+10])
        xdotvelocity_list.append(xdotvelocity[i:i+10])
        ydotvelocity_list.append(ydotvelocity[i:i+10])
        zdotvelocity_list.append(zdotvelocity[i:i+10])
        corefactor_list.append(corefactor[i:i+10])
        sigmaFCN_list.append(sigmaFCN[i:i+10])
        xLaRa_list.append(xLaRa[i:i+10])
        yLaRa_list.append(yLaRa[i:i+10])
        zLaRa_list.append(zLaRa[i:i+10])
        cos1spin_list.append(cos1spin[i:i+10])
        sin1spin_list.append(sin1spin[i:i+10])
        cos2spin_list.append(cos2spin[i:i+10])
        sin2spin_list.append(sin2spin[i:i+10])
        cos3spin_list.append(cos3spin[i:i+10])
        sin3spin_list.append(sin3spin[i:i+10])
        cos4spin_list.append(cos4spin[i:i+10])
        sin4spin_list.append(sin4spin[i:i+10])
        xpcos1_list.append(xpcos1[i:i+10])
        xpsin1_list.append(xpsin1[i:i+10])
        ypcos1_list.append(ypcos1[i:i+10])
        ypsin1_list.append(ypsin1[i:i+10])
        xpcos2_list.append(xpcos2[i:i+10])
        xpsin2_list.append(xpsin2[i:i+10])
        ypcos2_list.append(ypcos2[i:i+10])
        ypsin2_list.append(ypsin2[i:i+10])
        xpcos3_list.append(xpcos3[i:i+10])
        xpsin3_list.append(xpsin3[i:i+10])
        ypcos3_list.append(ypcos3[i:i+10])
        ypsin3_list.append(ypsin3[i:i+10])
        xpcos4_list.append(xpcos4[i:i+10])
        xpsin4_list.append(xpsin4[i:i+10])
        ypcos4_list.append(ypcos4[i:i+10])
        ypsin4_list.append(ypsin4[i:i+10])
        xpcos5_list.append(xpcos5[i:i+10])
        xpsin5_list.append(xpsin5[i:i+10])
        ypcos5_list.append(ypcos5[i:i+10])
        ypsin5_list.append(ypsin5[i:i+10])
        i = i+10

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

    # yposition plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],yposition_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $y$ [m]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/yposition_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # zposition plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],zposition_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $z$ [m]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/zposition_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # xdotvelocity plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],xdotvelocity_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $\dot{x}$ [m/s]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xdotvelocity_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # ydotvelocity plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],ydotvelocity_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $\dot{y}$ [m/s]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ydotvelocity_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # zdotvelocity plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],zdotvelocity_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $\dot{z}$ [m/s]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/zdotvelocity_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # corefactor plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],corefactor_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $F$ [-]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/corefactor_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # sigmaFCN plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],sigmaFCN_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $\sigma_{FCN}$ [rad/s]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/sigmaFCN_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # xLaRa plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],xLaRa_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $x_{LaRa}$ [m]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xLaRa_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # yLaRa plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],yLaRa_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $y_{LaRa}$ [m]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/yLaRa_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # zLaRa plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],zLaRa_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $z_{LaRa}$ [m]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/zLaRa_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    plt.rcParams.update({'font.size': 20})
    # cos1spin plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],cos1spin_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $\phi^c_1$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/cos1spin_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # sin1spin plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],sin1spin_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $\phi^s_1$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/sin1spin_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # cos2spin plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],cos2spin_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $\phi^c_2$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/cos2spin_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # sin2spin plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],sin2spin_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $\phi^s_2$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/sin2spin_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # cos3spin plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],cos3spin_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $\phi^c_3$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/cos3spin_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # sin3spin plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],sin3spin_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $\phi^s_3$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/sin3spin_plot.pdf"),bbox_inches="tight")
    plt.show() 
    plt.close("all") 

    # cos4spin plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],cos4spin_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $\phi^c_4$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/cos4spin_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # sin4spin plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],sin4spin_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ $\phi^s_4$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/sin4spin_plot.pdf"),bbox_inches="tight")
    plt.show() 
    plt.close("all")

    # xpcos1 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],xpcos1_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${X_p}^c_1$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpcos1_plot.pdf"),bbox_inches="tight")
    plt.show() 
    plt.close("all")

    # xpsin1 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],xpsin1_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${X_p}^s_1$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpsin1_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # ypcos1 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],ypcos1_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${Y_p}^c_1$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypcos1_plot.pdf"),bbox_inches="tight")
    plt.show() 
    plt.close("all")

    # ypsin1 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],ypsin1_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${Y_p}^s_1$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypsin1_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # xpcos2 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],xpcos2_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${X_p}^c_2$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpcos2_plot.pdf"),bbox_inches="tight")
    plt.show() 
    plt.close("all")

    # xpsin2 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],xpsin2_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${X_p}^s_2$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpsin2_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # ypcos2 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],ypcos2_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${Y_p}^c_2$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypcos2_plot.pdf"),bbox_inches="tight")
    plt.show() 
    plt.close("all")

    # ypsin2 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],ypsin2_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${Y_p}^s_2$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypsin2_plot.pdf"),bbox_inches="tight")
    plt.show() 
    plt.close("all")

    # xpcos3 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],xpcos3_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${X_p}^c_3$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpcos3_plot.pdf"),bbox_inches="tight")
    plt.show() 
    plt.close("all")

    # xpsin3 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],xpsin3_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${X_p}^s_3$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpsin3_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # ypcos3 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],ypcos3_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${Y_p}^c_3$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypcos3_plot.pdf"),bbox_inches="tight")
    plt.show() 
    plt.close("all")

    # ypsin3 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],ypsin3_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${Y_p}^s_3$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypsin3_plot.pdf"),bbox_inches="tight")
    plt.show() 
    plt.close("all")

    # xpcos4 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],xpcos4_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${X_p}^c_C$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpcos4_plot.pdf"),bbox_inches="tight")
    plt.show() 
    plt.close("all")

    # xpsin4 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],xpsin4_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${X_p}^s_C$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpsin4_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # ypcos4 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],ypcos4_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${Y_p}^c_C$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypcos4_plot.pdf"),bbox_inches="tight")
    plt.show() 
    plt.close("all")

    # ypsin4 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],ypsin4_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${Y_p}^s_C$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypsin4_plot.pdf"),bbox_inches="tight")
    plt.show() 
    plt.close("all")

    # xpcos5 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],xpcos5_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${X_p}^c_4$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpcos5_plot.pdf"),bbox_inches="tight")
    plt.show() 
    plt.close("all")

    # xpsin5 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],xpsin5_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${X_p}^s_4$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/xpsin5_plot.pdf"),bbox_inches="tight")
    plt.show()
    plt.close("all")

    # ypcos5 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],ypcos5_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${Y_p}^c_4$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypcos5_plot.pdf"),bbox_inches="tight")
    plt.show() 
    plt.close("all")

    # ypsin5 plot
    plt.figure(figsize=(6,6))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    plt.gca().ticklabel_format(useOffset=False)
    for i in range(0,10):
        plt.plot(corr_list[i],ypsin5_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    #plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    plt.ylabel(r"1-$\sigma$ ${Y_p}^s_4$ [mas]")
    plt.grid()
    plt.colorbar(mpl.cm.ScalarMappable(norm=mpl.colors.Normalize(vmin=1, vmax=10),cmap=mpl.cm.jet),label='Number of Stations')
    plt.savefig(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/vandv/ypsin5_plot.pdf"),bbox_inches="tight")
    plt.show() 
    plt.close("all")

print("--- %s seconds ---" % (time.time() - run_time))