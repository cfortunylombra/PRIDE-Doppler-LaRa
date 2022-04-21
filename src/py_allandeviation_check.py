"""
Description: Allan Deviation Check

Author: C. Fortuny-Lombraña
"""

import time

run_time = time.time()
if __name__=="__main__":
    ########################################################################################################################
    ################################################## IMPORT PACKAGES #####################################################
    ########################################################################################################################

    import os
    import json
    import copy
    import matplotlib.pyplot as plt
    import allantools
    import numpy as np
    import numpy as np
    from scipy import stats
    from scipy.stats import norm
    from collections import Counter

    import sys
    sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")
    from tudatpy.kernel import constants

    ########################################################################################################################
    ################################################## IMPORT OUTPUT FILES #################################################
    ########################################################################################################################

    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/POD_LaRa')
    os.makedirs(output_folder_path,exist_ok=True)

    # Import the first observation (time list)
    input = np.loadtxt(output_folder_path+'/concatenated_times.dat',delimiter=' ')
    time_array = input[:50]

    # Import the first observation (doppler residual list)
    input = np.loadtxt(output_folder_path+'/doppler_residuals.dat',delimiter=' ')
    doppler_last_iteration = input[:50,-1]*constants.SPEED_OF_LIGHT_LONG
    
    # Plot the Allan deviation using real-data
    plt.figure(figsize=(7,7))
    rate = 1/Counter(np.diff(time_array)).most_common(1)[0][0]
    print(rate)
    taus2, adevs, errors, ns = allantools.mdev(np.array(doppler_last_iteration),\
        rate = rate, data_type = 'freq',taus='decade')
    adevs_white = list()
    adevs_white.append(adevs[0]-errors[0])
    for i in range(1,len(taus2)):
        adevs_white.append(10**(-1/2*(np.log10(taus2[i])-np.log10(taus2[i-1]))+np.log10(adevs_white[i-1])))
    plot1 = plt.errorbar(taus2,adevs,yerr=errors,ecolor='green')
    plot2 = plt.errorbar(taus2, adevs_white,ecolor='red')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Tau [s]')
    plt.ylabel('Allan Deviation')
    plt.legend([plot1,plot2],['Real-data','Simulated White Noise'])
    plt.grid()
    plt.axis('equal')
    plt.show()
    plt.close('all')


print("--- %s seconds ---" % (time.time() - run_time))