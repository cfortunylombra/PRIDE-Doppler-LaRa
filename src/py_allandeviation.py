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
    import glob
    import numpy as np
    #Install pip install pandas
    import pandas as pd
    import math
    import datetime
    import json
    import matplotlib.pyplot as plt

    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/InSight')
    os.makedirs(output_folder_path,exist_ok=True)

    with open(output_folder_path+'/Fdets_data.json', 'r') as fp:
        data_fdets = json.load(fp)

    with open(output_folder_path+'/Phases_data.json', 'r') as fp:
        data_phases = json.load(fp)

    for fdets_station_pointer in data_fdets.keys():
        for fdets_station_starttime_pointer in data_fdets[fdets_station_pointer].keys():
            plt.plot(data_fdets[fdets_station_pointer][fdets_station_starttime_pointer]['Time(UTC) [s]'],\
                data_fdets[fdets_station_pointer][fdets_station_starttime_pointer]['Doppler noise [Hz]'])
            plt.show()

print("--- %s seconds ---" % (time.time() - run_time))