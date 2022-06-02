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

    np.set_printoptions(suppress=False,precision=10)
    corefactor = np.loadtxt(os.path.dirname(os.path.realpath(__file__)).replace('/src',"/output/verification_cos4spin.dat"))

    corr_list = list()
    station_number_list = list()
    corefactor_list = list()

    i = 0
    for station in np.linspace(1,10,10):
        corr_list.append(np.linspace(0,9,10)/10)
        station_number_list.append(station*np.ones(10))
        corefactor_list.append(corefactor[i:i+10])
        i = i+10

    plt.figure(figsize=(10,10))
    colormap = plt.cm.gist_ncar
    plt.gca().set_prop_cycle(plt.cycler('color', plt.cm.jet(np.linspace(0, 1, 10))))
    for i in range(0,10):
        plt.plot(corr_list[i],corefactor_list[i],'-o',label=str(int(station_number_list[i][0]))+" Station/s")
    plt.legend()
    plt.xlabel(r"Correlation $\rho$ [-]")
    #plt.ylabel(r"1-$\sigma$ F [-]")
    plt.ylabel(r"1-$\Phi^c_4$ [mas]")
    plt.grid()
    plt.show()
    

print("--- %s seconds ---" % (time.time() - run_time))