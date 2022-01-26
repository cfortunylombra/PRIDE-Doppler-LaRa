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
    #Install pip install pandas
    import json
    import matplotlib.pyplot as plt

    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/InSight')
    os.makedirs(output_folder_path,exist_ok=True)

    ########################################################################################################################
    ################################################## PLOTS ###############################################################
    ########################################################################################################################

    # Open the Fdets JSON file
    with open(output_folder_path+'/Fdets_data.json', 'r') as fp:
        data_fdets = json.load(fp)

    # Choose x- and y-axis (Comment out the useful x_axis and y_axis)
    x_axis_fdets = 'Time(UTC) [s]'
    #y_axis_fdets = 'Doppler noise [Hz]'
    #y_axis_fdets = 'Signal-to-Noise ratio'
    y_axis_fdets = 'Spectral max'
    #y_axis_fdets = 'Freq. detection [Hz]'

    boolean_fdets = False
    if boolean_fdets:
        # Iterate along the ground stations inside the Fdets JSON file
        for fdets_station_pointer in data_fdets.keys():
            # Iterate along the first time stamp for each ground station
            for fdets_station_starttime_pointer in data_fdets[fdets_station_pointer].keys():
                plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][x_axis_fdets],\
                    data_fdets[fdets_station_pointer][fdets_station_starttime_pointer][y_axis_fdets])
                plt.xlabel(x_axis_fdets)
                plt.ylabel(y_axis_fdets)
                plt.title(fdets_station_pointer+' station at '+fdets_station_starttime_pointer)
                plt.show()

    # Open the Phases JSON file
    with open(output_folder_path+'/Phases_data.json', 'r') as fp:
        data_phases = json.load(fp)

    # Choose x- and y-axis (Comment out the useful x_axis and y_axis)
    x_axis_phases = 'Time stamp [s]'
    y_axis_phases = 'Phase [rad]'

    boolean_phases = True
    if boolean_phases:
        # Iterate along the ground stations inside the Phases JSON file
        for phases_station_pointer in data_phases.keys():
            # Iterate along the first time stamp for each ground station
            for phases_station_starttime_pointer in data_phases[phases_station_pointer].keys():
                plt.plot(data_phases[phases_station_pointer][phases_station_starttime_pointer][x_axis_phases],\
                    data_phases[phases_station_pointer][phases_station_starttime_pointer][y_axis_phases])
                plt.xlabel(x_axis_phases)
                plt.ylabel(y_axis_phases)
                plt.title(phases_station_pointer+' station at '+phases_station_starttime_pointer)
                plt.show()

print("--- %s seconds ---" % (time.time() - run_time))