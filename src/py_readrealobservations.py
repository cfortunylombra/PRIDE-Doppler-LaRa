"""
Description: Read Real Observations

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
    import pandas as pd # This package might be not installed, please `pip install pandas`
    import math
    import datetime
    import json

    ########################################################################################################################
    ################################################## FUNCTIONS ###########################################################
    ########################################################################################################################

    def mjd_to_jd(mjd):
        """
        Convert Modified Julian Day to Julian Day.
        Inputs:
        mjd : float
            Modified Julian Day
        Returns:
        jd : float
            Julian Day   
        """
        return mjd + 2400000.5

    def jd_to_date(jd):
        """
        Convert Julian Day to date.
        Algorithm from 'Practical Astronomy with your Calculator or Spreadsheet', 
            4th ed., Duffet-Smith and Zwart, 2011.
        Inputs:
        jd : float
            Julian Day   
        Returns:
        year : int
            Year as integer. Years preceding 1 A.D. should be 0 or negative.
            The year before 1 A.D. is 0, 10 B.C. is year -9.
        month : int
            Month as integer, Jan = 1, Feb. = 2, etc.
        day : float
            Day, may contain fractional part.       
        """
        jd = jd + 0.5
        F, I = math.modf(jd)
        I = int(I)
        A = math.trunc((I - 1867216.25)/36524.25)
        if I > 2299160:
            B = I + 1 + A - math.trunc(A / 4.)
        else:
            B = I
        C = B + 1524
        D = math.trunc((C - 122.1) / 365.25)
        E = math.trunc(365.25 * D)
        G = math.trunc((C - E) / 30.6001)
        day = C - E + F - math.trunc(30.6001 * G)
        if G < 13.5:
            month = G - 1
        else:
            month = G - 13
        if month > 2.5:
            year = D - 4716
        else:
            year = D - 4715
        return year, month, int(day)

    def UTC_to_h_min_s(utc):
        """
        Convert UTC [s] to Time in hours. minutes and seconds.
        Inputs:
        utc : float
            UTC [s]
        Returns:
        hour : float
             hours
        min : float
             minutes
        sec : float
             seconds  
        """
        hour = utc//3600%24
        min = utc//60%60
        sec = utc%60
        return (int(hour), int(min), int(sec))

    ########################################################################################################################
    ################################################## READING FILES #######################################################
    ########################################################################################################################

    # Initialize dataset
    data_fdets = dict()
    data_phases = dict()

    # List all the folders inside the ED045 folder
    folders_per_time_scan = glob.glob(os.path.dirname(os.path.realpath(__file__))+'/ED045/*')

    # Iterate along the folders inside the ED045 folder
    for folder_time_scan_pointer in range(0,len(folders_per_time_scan)):

    ########################################################################################################################
    ################################################## FDETS FILES #######################################################
    ########################################################################################################################

        # List all the files that start with Fdets 
        files_fdets = glob.glob(folders_per_time_scan[folder_time_scan_pointer]+'/Fdets*.txt')

        # Iterate along the Fdets files 
        for fdets_pointer in range(0,len(files_fdets)):
            #print(files_fdets[fdets_pointer])
            file_name = os.path.basename(files_fdets[fdets_pointer])
            file_name_split = file_name.split('.')
            station_name = file_name_split[4]
            
            # Read each file
            with open(files_fdets[fdets_pointer]) as file:
                lines = file.read().splitlines()

                # Count the skip rows with counting the times of the # element
                skip_rows = 0

                # Read each line
                for pointer_line in range(0,len(lines)):
                    line = lines[pointer_line]
                    # Split the line
                    line_split = line.split()

                    # Counting how many skip rows
                    if line.count('#') == 1:
                        skip_rows +=1
                    
                    # Taking the information about base frequency, dF and dT
                    if pointer_line == 1:
                        base_frequency_mhz = float(line_split[line_split.index('frequency:')+1:line_split.index('frequency:')+2][0])
                        #print(base_frequency_mhz)
                        boolean_df = False
                        boolean_dt = False
                        if line_split.count('dF:') == 1:
                            df_hz = float(line_split[line_split.index('dF:')+1:line_split.index('dF:')+2][0])
                            #print(df_hz)
                            boolean_df = True
                        if line_split.count('dT:') == 1:
                            dt_s = float(line_split[line_split.index('dT:')+1:line_split.index('dT:')+2][0])
                            #print(dt_s)
                            boolean_dt = True

                    # Taking the column names as a list for the dictionary
                    if pointer_line == 2:
                        if ':' in line_split:
                            columns = line_split[line_split.index(':')+1:]
                        elif 'Format:' in line_split:
                            columns = line_split[line_split.index('Format:')+1:]
                        column_names = list()
                        column_element = ''
                        for string_pointer in range(0,len(columns)):
                            if columns[string_pointer] == '|':
                                column_names.append(column_element)
                                column_element = ''
                            elif string_pointer == len(columns)-1:
                                column_element +=  ' ' + columns[string_pointer]
                                column_names.append(column_element)
                            else:
                                if column_element == '':
                                    column_element = columns[string_pointer]
                                else:
                                    column_element +=  ' ' + columns[string_pointer]

            # Build a pandas frame for each file
            dataframe_fdets = pd.read_csv(files_fdets[fdets_pointer], sep=' ', names=column_names,skiprows = skip_rows)

            # Date tuple (year, month, day)
            date_tuple = jd_to_date(mjd_to_jd(dataframe_fdets[column_names[0]][0]))

            # Time tuple (hours, minutes, seconds)
            time_tuple = UTC_to_h_min_s(dataframe_fdets[column_names[1]][0])

            # Create a time stamp
            datetime_tuple = date_tuple+time_tuple
            time_stamp = datetime.datetime(datetime_tuple[0],datetime_tuple[1],datetime_tuple[2],datetime_tuple[3],datetime_tuple[4],datetime_tuple[5])
            #print(time_stamp)

            # Insert all the information into a dictionary
            if not station_name in data_fdets.keys():
                data_fdets[station_name] = dict()
            data_fdets[station_name][str(time_stamp)] = dict()

            for column_pointer in range(0,len(column_names)):
                if column_names[column_pointer] == 'Signal-to-Noise ratio':
                    data_fdets[station_name][str(time_stamp)][column_names[column_pointer]+' [dB]']=list(10*np.log10(dataframe_fdets[column_names[column_pointer]]))
                else:
                    data_fdets[station_name][str(time_stamp)][column_names[column_pointer]]=list(dataframe_fdets[column_names[column_pointer]])
            data_fdets[station_name][str(time_stamp)]['Base Frequency [MHz]']=base_frequency_mhz
            if boolean_df:
                data_fdets[station_name][str(time_stamp)]['dF [Hz]']=df_hz
            if boolean_dt:
                data_fdets[station_name][str(time_stamp)]['dT [s]']=dt_s

    ########################################################################################################################
    ################################################## PHASES FILES #######################################################
    ########################################################################################################################

        # Iterate along the Phases files 
        files_phases = glob.glob(folders_per_time_scan[folder_time_scan_pointer]+'/Phases*.txt')
        for phases_pointer in range(0,len(files_phases)):
            #print(files_phases[phases_pointer])
            file_name = os.path.basename(files_phases[phases_pointer])
            file_name_split = file_name.split('.')
            station_name = file_name_split[4]

            # Read each file
            with open(files_phases[phases_pointer]) as file:
                lines = file.read().splitlines()

                # Count the skip rows with counting the times of the # element
                skip_rows = 0

                # Read each line
                for pointer_line in range(0,len(lines)):
                    line = lines[pointer_line]
                    # Split the line
                    line_split = line.split()

                    # Counting how many skip rows
                    if line.count('#') == 1:
                        skip_rows +=1

                    # Taking the column names as a list for the dictionary
                    if pointer_line == 1:
                        labels = line_split[line_split.index('#')+1:]
                        labels_names = list()
                        label_element = ''
                        for string_pointer in range(0,len(labels)):
                            if labels[string_pointer] == '-':
                                labels_names.append(label_element)
                                label_element = ''
                            elif string_pointer == len(labels)-1:
                                label_element +=  ' ' + labels[string_pointer]
                                labels_names.append(label_element)
                            else:
                                if label_element == '':
                                    label_element = labels[string_pointer]
                                else:
                                    label_element +=  ' ' + labels[string_pointer]

                    # Taking the values as a list for the dictionary
                    if pointer_line == 2:
                        specs = line_split[line_split.index('#')+1:]
                        specs[0] = specs[0].replace('\\', '')
                        specs = list(np.float_(specs))

                    # Taking the column names as a list for the dictionary
                    if pointer_line == 3:
                        columns = line_split[line_split.index('#')+1:]
                        column_names = list()
                        column_element = ''
                        for string_pointer in range(0,len(columns)):
                            if columns[string_pointer] == '|':
                                column_names.append(column_element)
                                column_element = ''
                            elif string_pointer == len(columns)-1:
                                column_element +=  ' ' + columns[string_pointer]
                                column_names.append(column_element)
                            else:
                                if column_element == '':
                                    column_element = columns[string_pointer]
                                else:
                                    column_element +=  ' ' + columns[string_pointer]

            # Build a pandas frame for each file
            dataframe_phases = pd.read_csv(files_phases[phases_pointer], sep=' ', names=column_names,skiprows = skip_rows)
            dataframe_phases[column_names[0]] +=specs[1]
             
            # Date tuple (year, month, day)
            date_tuple = jd_to_date(mjd_to_jd(specs[0]))

            # Time tuple (hours, minutes, seconds)
            time_tuple = UTC_to_h_min_s(specs[1])

            # Create a time stamp
            datetime_tuple = date_tuple+time_tuple
            time_stamp = datetime.datetime(datetime_tuple[0],datetime_tuple[1],datetime_tuple[2],datetime_tuple[3],datetime_tuple[4],datetime_tuple[5])          
            #print(time_stamp)

            # Insert all the information into a dictionary
            if not station_name in data_phases.keys():
                data_phases[station_name] = dict()
            data_phases[station_name][str(time_stamp)] = dict()

            for column_pointer in range(0,len(column_names)):
                data_phases[station_name][str(time_stamp)][column_names[column_pointer]]=list(dataframe_phases[column_names[column_pointer]])
            
            for specs_pointer in range(0,len(specs)):
                data_phases[station_name][str(time_stamp)][labels_names[specs_pointer]]=specs[specs_pointer]
                
    # Save all the data into JSON files
    output_folder_path = os.path.dirname(os.path.realpath(__file__)).replace('/src','/output/InSight')
    os.makedirs(output_folder_path,exist_ok=True)

    with open(output_folder_path+'/Fdets_data.json', 'w') as fp:
        json.dump(data_fdets, fp)

    with open(output_folder_path+'/Phases_data.json', 'w') as fp:
        json.dump(data_phases, fp)

print("--- %s seconds ---" % (time.time() - run_time))