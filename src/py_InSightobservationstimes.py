import time

run_time = time.time()
if __name__=="__main__":

    import os
    import numpy as np
    import datetime

    transmitter_ground_station_name = "DSS 63" # DSS 14

    with open(os.path.dirname(os.path.realpath(__file__))+'/nsyt_maro_2018_331_2020_366.tdm') as file:
        lines = file.read().splitlines()
        meta_start_line = lines.index('PARTICIPANT_1          = '+transmitter_ground_station_name)
        data_start_line = list(filter(lambda x: x >= meta_start_line, [i for i, x in enumerate(lines) if x == 'DATA_START']))[0]+1
        data_stop_line = list(filter(lambda x: x >= meta_start_line, [i for i, x in enumerate(lines) if x == 'DATA_STOP']))[0] 
        lines = lines[data_start_line:data_stop_line]

        datetime_list = list()
        epoch_list = list()
        for pointer_time in range(0,len(lines)):
            line = lines[pointer_time]
            date_and_time = line.split()[2]
            year = int(date_and_time.split('-')[0])
            day_of_year = int((date_and_time.split('-')[1]).split('T')[0])
            hour_min_sec = ((date_and_time.split('-')[1]).split('T')[1]).split(':')
            hour = int(hour_min_sec[0])
            min = int(hour_min_sec[1])
            sec = int((hour_min_sec[2].split('.'))[0])

            date = datetime.datetime(year,1,1)+datetime.timedelta(days=day_of_year-1,hours=hour,minutes=min,seconds=sec)
            datetime_list.append(date)

            #epoch = (date - datetime.datetime(2000,1,1,12,0,0)).total_seconds()
            epoch = (date - datetime.datetime(2000,1,1,11,58,55,816000)).total_seconds()
            epoch_list.append(epoch)

    print(datetime_list[0],epoch_list[0])

    with open(os.path.dirname(os.path.realpath(__file__))+'/InSight_mes_upto_31122021.forCarlos') as file:
        lines = file.read().splitlines()

        transmitter_ground_station_number =  [int(s) for s in transmitter_ground_station_name.split() if s.isdigit()][0]
        
        datetime_list = list()
        epoch_list = list()
        for pointer_time in range(0,len(lines)):
            line = lines[pointer_time]
            line_info = line.split()
            if float(line_info[0]) == transmitter_ground_station_number and float(line_info[1]) == transmitter_ground_station_number:
                #epoch_list.append(float(line_info[2]))
                epoch_list.append(float(line_info[2])+datetime.timedelta(minutes=1,seconds=4,milliseconds=184).total_seconds())

                date = datetime.datetime(2000,1,1,12,0,0)+datetime.timedelta(seconds=float(line_info[2]))
                datetime_list.append(date)

    print(datetime_list[0],epoch_list[0])

print("--- %s seconds ---" % (time.time() - run_time))