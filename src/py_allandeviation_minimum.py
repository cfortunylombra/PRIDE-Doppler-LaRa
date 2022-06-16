"""
Description: Save the minimum Allan Deviation for InSight observation span

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
    from scipy import stats

    base_frequency = 8400.5*10**6 #Hz

    station_nickname = {"Mc":"MEDICINA","Wz":"WETTZELL","O6":"ONSALA60","Ef":"EFLSBERG","Wb":"WRT0","Ys":"YEBES40M",\
        "Ho":"HOBART26","T6":"TIANMA65","Bd":"BADARY","Ww":"WARK30M","Cd":"CEDUNA","Hh":"HARTRAO","Ir":"IRBENE"}
    
    minimum_ftes_allan_deviation = dict()

    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/InSight')
    os.makedirs(output_folder_path,exist_ok=True)

    # Open the Fdets JSON file
    with open(output_folder_path+'/Fdets_data_updated.json', 'r') as fp:
        data_fdets = json.load(fp)

    for fdets_station_pointer in data_fdets.keys():
        allan_deviations_per_station = []
        for fdets_station_starttime_pointer in data_fdets[fdets_station_pointer].keys():
            allan_deviations_per_station.append(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer]['Sigma'][0][0])
        if 0 in allan_deviations_per_station:
            allan_deviations_per_station.remove(0)
        minimum_ftes_allan_deviation[station_nickname[fdets_station_pointer]] = min(allan_deviations_per_station)

    minimum_ftes_allan_deviation["DSN"] = 0.21545336*10**(-3)/base_frequency

    print(minimum_ftes_allan_deviation)

    print("Common noise: ",2.08650e-14)

    maximum_SNR = dict()

    # Open the Fdets JSON file
    with open(output_folder_path+'/Fdets_data_updated.json', 'r') as fp:
        data_fdets = json.load(fp)

    for fdets_station_pointer in data_fdets.keys():
        SNR_per_station = []
        for fdets_station_starttime_pointer in data_fdets[fdets_station_pointer].keys():
            SNR_per_station.extend(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer]["Signal-to-Noise ratio [dB]"][0]) 
        maximum_SNR[station_nickname[fdets_station_pointer]] = max(SNR_per_station)

    print(maximum_SNR)

    minimum_phases_allan_deviation = dict()

    with open(output_folder_path+'/Phases_data_updated.json', 'r') as fp:
        data_phases = json.load(fp)

    for phases_station_pointer in data_phases.keys():
        allan_deviations_per_station = []
        for phases_station_starttime_pointer in data_phases[phases_station_pointer].keys():
            allan_deviations_per_station.append(data_phases[phases_station_pointer][phases_station_starttime_pointer]['Sigma'])
        if 0 in allan_deviations_per_station:
            allan_deviations_per_station.remove(0)
        minimum_phases_allan_deviation[station_nickname[phases_station_pointer]] = min(allan_deviations_per_station)

    print(minimum_phases_allan_deviation)

    minimum_allan_deviation = dict()

    for fdets_station_pointer in data_fdets.keys():
        if fdets_station_pointer in data_phases.keys():
            minimum_allan_deviation[station_nickname[fdets_station_pointer]] = min(minimum_ftes_allan_deviation[station_nickname[fdets_station_pointer]],
                minimum_phases_allan_deviation[station_nickname[fdets_station_pointer]])
        else:
            minimum_allan_deviation[station_nickname[fdets_station_pointer]] = minimum_ftes_allan_deviation[station_nickname[fdets_station_pointer]]

    print(minimum_allan_deviation)

print("--- %s seconds ---" % (time.time() - run_time))