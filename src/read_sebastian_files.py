"""
Description: Read Sebastian files

Author: C. Fortuny-Lombraña
"""

from selectors import EpollSelector
import time

import scipy

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
    import scipy
    import scipy.fftpack
    import pandas as pd

    import sys
    sys.path.insert(0, "/home/cfortunylombra/tudat-bundle/cmake-build-release-wsl/tudatpy/")
    from tudatpy.kernel import constants


    # Change directory
    output_folder_path = os.path.dirname(os.path.realpath(__file__))
    os.makedirs(output_folder_path,exist_ok=True)

    # Remove warnings
    np.warnings.filterwarnings('ignore')

    time_days = list()
    mean_mHz = list()
    std_mHz = list()

    # Take all the means, standard deviations and times and append them to the empty lists
    with open(output_folder_path+'/ResStatPerPass_ForCarlos.txt') as f:
        lines = f.readlines()
        index_gap = None
        for line in lines[1:]:
            line_split = line.split()
            if not (np.isnan(float(line_split[1])) and np.isnan(float(line_split[2]))):
                time_days.append(float(line_split[0])) #*constants.JULIAN_DAY
                mean_mHz.append(float(line_split[1]))
                std_mHz.append(float(line_split[2]))
            else:
                index_gap = len(time_days)
    print(np.mean(std_mHz))
    # Plot - mean as a function of time
    plt.figure(figsize=(15,6))
    plt.plot(time_days,mean_mHz)
    plt.ylabel('Mean [mHz]')
    plt.xlabel('Time [days]')
    plt.grid()
    plt.show()
    plt.close('all')

    # Plot - standard deviation as a function of time
    plt.figure(figsize=(15,6))
    plt.plot(time_days,std_mHz)
    plt.ylabel('Std [mHz]')
    plt.xlabel('Time [days]')
    plt.grid()
    plt.show()
    plt.close('all')

    # Plot - mean+standard deviation as a function of time
    plt.figure(figsize=(15,6))
    plt.errorbar(time_days,mean_mHz,yerr=std_mHz)
    plt.ylabel('Mean$\pm$Std [mHz]')
    plt.xlabel('Time [days]')
    plt.grid()
    plt.show()
    plt.close('all')


    # Insert white random noise in the gap and use cubic splines
    np.random.seed(42)

    poly = np.polyfit(time_days,std_mHz,8)
    poly_entire = np.poly1d(poly)(time_days)
    std_error = np.mean(np.abs(poly_entire-std_mHz))

    poly_gap = np.poly1d(poly)(np.arange(time_days[index_gap-1]+1,time_days[index_gap],1))
    poly_noise = np.random.normal(poly_gap,std_error)

    f0 = scipy.interpolate.CubicSpline(time_days[:index_gap],std_mHz[:index_gap])
    f1 = scipy.interpolate.CubicSpline(time_days[index_gap:],std_mHz[index_gap:])
    f2 = scipy.interpolate.CubicSpline(time_days,std_mHz)

    y_concatenate = np.concatenate((std_mHz[:index_gap],poly_noise,std_mHz[index_gap:]))
    time_concatenate = np.concatenate((time_days[:index_gap],np.arange(time_days[index_gap-1]+1,time_days[index_gap],1),time_days[index_gap:]))

    f_total = scipy.interpolate.CubicSpline(time_concatenate,y_concatenate)

    plt.figure(figsize=(15,6))
    plt.plot(time_days,std_mHz,'r')
    plt.plot(np.arange(time_days[0],time_days[-1],1),f_total(np.arange(time_days[0],time_days[-1],1)),'g')
    plt.plot(np.arange(time_days[0],time_days[-1],1),np.poly1d(poly)(np.arange(time_days[0],time_days[-1],1)),'k')
    plt.plot(np.arange(time_days[0],time_days[-1],1),f2(np.arange(time_days[0],time_days[-1],1)),'b')
    plt.ylabel('Std [mHz]')
    plt.xlabel('Time [days]')
    plt.grid()
    plt.show()
    plt.close('all')

    # Use the nearest interpolation

    f3 = scipy.interpolate.interp1d(time_days, std_mHz, fill_value='extrapolate', kind='nearest')

    plt.figure(figsize=(15,6))
    plt.plot(time_days,std_mHz,'r')
    plt.plot(np.arange(time_days[0],time_days[-1],1),f3(np.arange(time_days[0],time_days[-1],1)),'g')
    plt.ylabel('Std [mHz]')
    plt.xlabel('Time [days]')
    plt.grid()
    plt.show()
    plt.close('all')

print("--- %s seconds ---" % (time.time() - run_time))