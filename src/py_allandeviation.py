"""
Description: Allan Deviation

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
    import matplotlib.pyplot as plt
    import allantools
    import numpy as np

    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/InSight')
    os.makedirs(output_folder_path,exist_ok=True)

    ########################################################################################################################
    ################################################## PLOTS ###############################################################
    ########################################################################################################################

    # Open the Fdets JSON file
    with open(output_folder_path+'/Fdets_data.json', 'r') as fp:
        data_fdets = json.load(fp)

    # Choose x- and y-axis (Comment out the useful x_axis and y_axis), we can use as well 'Doppler noise [Hz]', 'Spectral max', 'Freq. detection [Hz]'
    x_axis_fdets = 'Time(UTC) [s]'
    y_axis_fdets = 'Signal-to-Noise ratio'
    x2_axis_fdets = 'Doppler noise [Hz]'

    boolean_fdets = False
    if boolean_fdets:
        # Iterate along the ground stations inside the Fdets JSON file
        for fdets_station_pointer in data_fdets.keys():
            # Iterate along the first time stamp for each ground station
            for fdets_station_starttime_pointer in data_fdets[fdets_station_pointer].keys():
                # First plot: y_axis_fdets vs x_axis_fdets
                plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets],\
                    data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets],'bo-')
                plt.xlabel(x_axis_fdets)
                plt.ylabel(y_axis_fdets)
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.show()

                # Second plot: y_axis_fdets vs x2_axis_fdets
                plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x2_axis_fdets],\
                    data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets],'rx')
                plt.xlabel(x2_axis_fdets)
                plt.ylabel(y_axis_fdets)
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.show()

                # Allan deviation using the x2_axis_fdets
                rate_fdets = 1/np.diff(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets])[0]
                taus2, adevs, errors, ns = allantools.adev(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x2_axis_fdets],\
                    rate = rate_fdets, data_type = 'freq')
                plt.loglog(taus2, adevs,'go-')
                plt.xlabel('Tau')
                plt.ylabel('Allan Deviation')
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.grid()
                plt.axis('equal')
                plt.show()
                    

    # Open the Phases JSON file
    with open(output_folder_path+'/Phases_data.json', 'r') as fp:
        data_phases = json.load(fp)

    # Choose x- and y-axis (Comment out the useful x_axis and y_axis)
    x_axis_phases = 'Time stamp [s]'
    y_axis_phases = 'Phase [rad]'

    boolean_phases = False
    if boolean_phases:
        # Iterate along the ground stations inside the Phases JSON file
        for phases_station_pointer in data_phases.keys():
            # Iterate along the first time stamp for each ground station
            for phases_station_starttime_pointer in data_phases[phases_station_pointer].keys():
                # First plot: y_axis_phases vs x_axis_phases 
                plt.plot(data_phases[phases_station_pointer][phases_station_starttime_pointer][x_axis_phases],\
                    data_phases[phases_station_pointer][phases_station_starttime_pointer][y_axis_phases],'b')
                plt.xlabel(x_axis_phases)
                plt.ylabel(y_axis_phases)
                plt.title(phases_station_pointer+' station at '+phases_station_starttime_pointer)
                plt.grid()
                plt.show()

                # Allan deviation using x_axis_phases
                rate_phases = 1/np.diff(data_phases[phases_station_pointer][phases_station_starttime_pointer][x_axis_phases])[0]
                freq = allantools.phase2frequency(data_phases[phases_station_pointer][phases_station_starttime_pointer][y_axis_phases],rate_phases)
                taus2, adevs, errors, ns = allantools.adev(freq,rate = rate_phases, data_type = 'freq')
                plt.loglog(taus2, adevs,'go-')
                plt.xlabel('Tau')
                plt.ylabel('Allan Deviation')
                plt.title(phases_station_pointer+' station at '+phases_station_starttime_pointer)
                plt.grid()
                plt.axis('equal')
                plt.show()

print("--- %s seconds ---" % (time.time() - run_time))