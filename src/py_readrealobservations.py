"""
Description: Read Real Observations

Author: C. Fortuny-Lombraña
"""

from ctypes import pointer
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

    # Initialize dataset
    data = dict()

    # List with all the folders inside the ED045 folder
    folders_per_time_scan = glob.glob(os.path.dirname(os.path.realpath(__file__))+'/ED045/*')

    # Iterate along the folders inside the ED045 folder
    for folder_time_scan_pointer in range(0,len(folders_per_time_scan)):
        files_fdets = glob.glob(folders_per_time_scan[folder_time_scan_pointer]+'/Fdets*.txt')

        # Iterate along the Fdets files 
        for fdets_pointer in range(0,len(files_fdets)):
            
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

                    # Taking the column names as a list for the dataframe
                    if skip_rows == 3:
                        columns = line_split[line_split.index(':')+1:]
                        column_names = list()
                        column_element = ""
                        for string_pointer in range(0,len(columns)):
                            if columns[string_pointer] == "|":
                                column_names.append(column_element)
                                column_element = ""
                            elif string_pointer == len(columns)-1:
                                column_element +=  ' ' + columns[string_pointer]
                                column_names.append(column_element)
                            else:
                                if column_element == "":
                                    column_element = columns[string_pointer]
                                else:
                                    column_element +=  ' ' + columns[string_pointer]

            # Build a pandas frame for each file
            dataframe = pd.read_csv(files_fdets[fdets_pointer], sep=" ", names=column_names,skiprows = skip_rows)
            print(dataframe)
            break
        break
     

print("--- %s seconds ---" % (time.time() - run_time))